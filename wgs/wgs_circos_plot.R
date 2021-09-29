##Circos plot function call
##Only autosomal chromosomes
circos.plot.fun <- function(circos.plot.list, chrX){
  if(chrX){
  circos.par(start.degree = 90, gap.after = c(rep(1, times = 22), 5))
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22, "X")), species = "hg38")
} else {
  circos.par(start.degree = 90, gap.after = c(rep(1, times = 21), 5))
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22)), species = "hg38")
  }
  # invisible(lapply(circos.plot.list, function(x) circos.genomicTrack(x, ylim=c(0,4), numeric.column = 4,
  #                                                             panel.fun = function(region, value, ...) {
  #                                                               i = getI(...)
  #                                                               circos.genomicPoints(region, value, 
  #                                                                                    col =ifelse(value >=1.2, "goldenrod",
  #                                                                                                ifelse(value <= 0.8, "royalblue", "gray")), 
  #                                                                                    pch = 16, cex = 0.3)
  #                                                               circos.lines(CELL_META$cell.xlim, c(1, 1), lty = 1, col = "gray")
  #                                                               for(h in seq(2, 4, by = 1)){
  #                                                                 circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "gray")}
  #                                                               circos.yaxis(side = "left", at = c(1, 2, 3, 4),sector.index = get.all.sector.index()[1], tick = T,
  #                                                                            labels.cex = 0.4, labels.niceFacing = T)
  #                                                             })))
  invisible(lapply(circos.plot.list, function(x) circos.genomicTrack(x, ylim=c(-2,5), numeric.column = 4,
                                                                     panel.fun = function(region, value, ...) {
                                                                       i = getI(...)
                                                                       circos.genomicPoints(region, value, 
                                                                                            col =ifelse(value >=0.32, "goldenrod",
                                                                                                        ifelse(value <= -0.41, "royalblue", "gray")), 
                                                                                            pch = 16, cex = 0.3)
                                                                       circos.lines(CELL_META$cell.xlim, c(0,0), lty = 1, col = "gray")
                                                                       for(h in seq(-2, 5, by = 2)){
                                                                         circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "gray")}
                                                                       circos.yaxis(side = "left", at = c(-2, 0, 1, 2, 3, 4,5),sector.index = get.all.sector.index()[1], tick = T,
                                                                                    labels.cex = 0.4, labels.niceFacing = T)
                                                                     })))
}