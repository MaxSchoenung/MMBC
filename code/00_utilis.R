add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

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


# Get TSS and transcripts -------------------------------------------------
tss_annot <- function(txdb,autosomes=T){
  tx_by_gn <- transcriptsBy(txdb, by="gene")
  unlisted <- unlist(tx_by_gn)
  TSS <- ifelse(strand(unlisted) == "+", start(unlisted), end(unlisted))
  TSS <- GRanges(seqnames(unlisted), IRanges(TSS, width=1), strand(unlisted))
  TSS$tx_name <- unlisted$tx_name
  if(autosomes==T){
    TSS <- TSS[seqnames(TSS)%in%paste0("chr",1:19),]
  }
  TSS <- sort(TSS)
  return(TSS)
}

# annotate_window(test.gr[1:12,],tss.mm10,5000)
annotate_window <- function(gr.query,tss,window=5000){
  list.ret <- list()
  for(i in 1:length(gr.query)){
    message(paste0("Processing ",i))
  overl.anno <- findOverlaps(gr.query[i,]+window,tss) 
  if(length(overl.anno)==0){
    list.ret[[i]] <-  NA
  }else{
    list.ret[[i]] <- tss[overl.anno@to,]$gene_id}
  }
  names(list.ret) <- names(gr.query)
  for(i in 1:length(list.ret)){
    list.ret[[i]] <- data.frame("site"=rep(names(list.ret)[i],length(list.ret[[i]])),
                                "gene"=list.ret[[i]],
                                "meth.diff"=rep(gr.query$mean.diff[i],length(list.ret[[i]])),
                                "CT"=rep(gr.query$CellType[i],length(list.ret[[i]]))
                                )
  }
  ret.df <- do.call(rbind,list.ret)
  return(ret.df)
}

