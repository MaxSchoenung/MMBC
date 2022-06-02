# Raster_Fun --------------------------------------------------------------
DrawGeomPointRast <- function(data, panel_params, coord, na.rm = FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300) {
  if (is.null(raster.width)) {
    raster.width <- par('fin')[1]
  }
  
  if (is.null(raster.height)) {
    raster.height <- par('fin')[2]
  }
  
  prev_dev_id <- dev.cur()
  
  p <- ggplot2::GeomPoint$draw_panel(data, panel_params, coord)
  dev_id <- Cairo::Cairo(type='raster', width=raster.width*raster.dpi, height=raster.height*raster.dpi, dpi=raster.dpi, units='px', bg="transparent")[1]
  
  grid::pushViewport(grid::viewport(width=1, height=1))
  grid::grid.points(x=p$x, y=p$y, pch = p$pch, size = p$size,name = p$name, gp = p$gp, vp = p$vp, draw = T)
  grid::popViewport()
  cap <- grid::grid.cap()
  dev.off(dev_id)
  dev.set(prev_dev_id)
  
  grid::rasterGrob(cap, x=0, y=0, width = 1,height = 1, default.units = "native",just = c("left","bottom"))
}

GeomPointRast <- ggplot2::ggproto(
  "GeomPointRast",
  ggplot2::GeomPoint,
  draw_panel = DrawGeomPointRast
)
