#dplyr, parallel, scales, tibble, hypertidy/vapour, hypertidy/ximage

# export VSI_CACHE=TRUE
# export GDAL_CACHEMAX=30%
# export VSI_CACHE_SIZE=10000000
# export GDAL_HTTP_MULTIPLEX=YES
# export GDAL_INGESTED_BYTES_AT_OPEN=32000
# export GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR
# export GDAL_HTTP_VERSION=2
# export GDAL_HTTP_MERGE_CONSECUTIVE_RANGES=YES
# export GDAL_NUM_THREADS=ALL_CPUS
#
system.time({
Sys.setenv(VSI_CACHE="TRUE")
Sys.setenv(GDAL_CACHEMAX="30%")
Sys.setenv(VSI_CACHE_SIZE="10000000")
Sys.setenv(GDAL_HTTP_MULTIPLEX="YES")
Sys.setenv(GDAL_INGESTED_BYTES_AT_OPEN="32000")
Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR")
Sys.setenv(GDAL_HTTP_VERSION="2")
Sys.setenv(GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES")
Sys.setenv(GDAL_NUM_THREADS="ALL_CPUS")

##' internal function to read a stac query and get all its hrefs (every variable)
hrefs0 <- function(x) {
  a <- jsonlite::fromJSON(x)
  l <- lapply(a$features$assets, \(.x) .x$href)
  nms <- names(a$features$assets)

  hrefs <- tibble::as_tibble(setNames(l, nms))
  if(("next" %in% a$links$rel)) {
    idx <- which("next" == a$links$rel)
    hrefs <- rbind(hrefs, Recall(a$links$href[idx]))
  }
  hrefs
}

## vectorized hrefs from a stac query
hrefs <- function(x) {
  out <- NULL
  for (x0 in x) out <- rbind(out, hrefs0(x0))
  out
}

## function to take a single-row of href, this works on red,green,blue,cloud from sentinel-2-c1-l2a
## returns a table of cell, red, green, blue - we mask cloud to 10 atm
warpfun <- function(x,  dim = c(1280, 0), ext = NULL, crs = NULL) {
  if (is.null(crs)) crs <- ""
  cl_arg <- NULL
  if (!is.null(dim)) cl_arg <- c(cl_arg, "-ts", dim[1], dim[2])
  if (!is.null(ext)) cl_arg <- c(cl_arg, "-te", ext[1], ext[3], ext[2], ext[4])
## I don't think we have to deal with bands, we get all or 1

  gdalraster::warp(x, tf <- tempfile(fileext = ".tif", tmpdir = "/vsimem"), t_srs = crs, cl_arg = cl_arg, quiet = TRUE)

  ds <- new(gdalraster::GDALRaster, tf)
  dat <- gdalraster::read_ds(ds)
  ds$close()
  gdalraster::vsi_unlink(tf)
  dat

}


## function to run the warp function and trim out the result, we only store valid pixels
## cell is relevant to dim+ext+crs, and used as a grouping for taking median
sclfun <- function(src,  dim = c(1280, 0), ext = NULL, crs = NULL) {

  check <- gdalraster::buildVRT(visual <- tempfile(fileext = ".vrt", tmpdir = "/vsimem"),
                                  sprintf("/vsicurl/%s", c(src$red, src$green, src$blue, src$cloud)), cl_arg = "-separate", quiet = TRUE)


  vis <- warpfun(visual, dim = dim, ext = ext, crs = crs)
  mm <- matrix(vis, ncol = 4)
  bad <- mm[,4] > 10
  out <- mm[,1:3]
  if (any(bad)) out[bad,] <- NA
  keep <- rowSums(out, na.rm = TRUE) > 0
  tibble::tibble(cell = which(keep),
                 red = out[keep,1], green = out[keep, 2], blue = out[keep, 3])

}


ql <- function(x, dim = NULL) {
  if (!grepl("/vsicurl", x)) x <- sprintf("/vsicurl/%s", x)
  info <- vapour::vapour_raster_info(x)
  dim <- tail(info$overviews, 2)
  ximage::ximage(d <- gdal_raster_data(x[1], target_dim = dim, bands = 1:info$bands), asp = 1)
  invisible(d)
}


## generate the scenes for 2023 for Fiji (we get two for each side of the anti-meridian)
#sds::stacit(c(176.1, 182.7, -20, -15), "2023")

ex <- c(-64, -63.5, -9, -8.5) ## we might use a different grid for the output, but this is the query

## this is Ryan's example comparing stackstac to odc-stac
## collection is sentinel-2-c1-l2a, if we use "sentinel-2-l2a" we have to use SCL not cloud
qu <- sds::stacit(ex, date = "2020-01")

## get the hrefs (26 for Ryan's example, I think c1 has cleaned up overlap a bit)
srcs <- hrefs(qu)

## specifying our output grid
ext <- ex
dm <- c(1280, 0)
crs <- "EPSG:4326"

## we need a template, these might be different so we normalize upfront
ref <- warpfun(srcs$aot[1], ext = ext, dim = dm, crs = crs)
ext <- attr(ref, "gis")$bbox[c(1, 3, 2, 4)]; dm <- attr(ref, "gis")$dim[1:2]; crs <- attr(ref, "gis")$srs

library(vapour)
library(ximage)
library(parallel)
n <- nrow(srcs)
cl <- makeCluster(n.cpus <- min(c(64, n)))
clusterExport(cl, "warpfun")
system.time(d <- parLapply(cl, split(srcs, 1:nrow(srcs))[sample(n)], sclfun, ext = ext, crs = crs))

stopCluster(cl)


library(dplyr)

## this is the slow part, takes 30 seconds
library(multidplyr)
system.time({
cluster <- new_cluster(18)
d <- bind_rows(d)
res0 <- d %>%   group_by(cell) %>% partition(cluster) %>% summarize(red = median(red), green = median(green),
                                                                               blue = median(blue), n = dplyr::n()) %>% collect()
rm(cluster)
})


dv <- range(unlist(res0[c("red", "green", "blue")]), na.rm = T)

res <- res0
clamp <- function(x, r) {
 x[x < r[1]] <- r[1]
 x[x > r[2]] <- r[2]
 x
}
res[c("red", "green", "blue")] <- lapply(res[c("red", "green", "blue")], \(.x) clamp(scales::rescale(.x, c(0, 255), from = c(350, 5000)), c(0, 255)))



atts <- attributes(ref)
ref <- replicate(3, rep(0.0, length(ref[[1]])), simplify = F)
attributes(ref) <- atts
ref[[1]][res$cell] <- as.integer(res$red)
ref[[2]][res$cell] <- as.integer(res$green)
ref[[3]][res$cell] <- as.integer(res$blue)

#png("sentinel.png", width = attr(ref, "dimension")[1], height = attr(ref, "dimension")[2])
par(mar = rep(0, 4))
ximage(ref, asp = 1/cos(mean(ext[3:4]) * pi/180), axes = FALSE, xlab = "", ylab = "")




refmask <- ref
refmask[[2]] <- refmask[[3]] <- NULL
refmask[[1]][res$cell] <- res$n
ximage(refmask, col = hcl.colors(12))

})
