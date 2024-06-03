#### batch functions 
Scatter_Density <- function(data = data, batch = batch, trt = NULL, expl.var = expl.var,
                            xlim = xlim, ylim = ylim, batch.legend.title = 'Batch', 
                            trt.legend.title = 'Treatment', density.lwd = 0.2,
                            title = NULL, title.cex = 1.5, legend.cex = 0.7, legend.title.cex =0.75){
  data = as.data.frame(data)
  batch = as.factor(batch)
  trt = as.factor(trt)
  if(nlevels(trt) >= 2){
    pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch, shape = trt)) + 
      geom_point() + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
      ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
      scale_color_manual(values = color.mixo(1:10)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title, shape = trt.legend.title)
    
    pTop <- ggplot(data,aes(x = data[ ,1], fill = batch, linetype = trt)) + 
      geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') + 
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
            plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
            axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(xlim[1], xlim[2]) + labs(title = title)
    
    pRight <- ggplot(data, aes(x=data[ ,2], fill = batch, linetype = trt)) + 
      geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
      theme(axis.title.x = element_text(size = rel(0.8)), 
            axis.title.y = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(ylim[1], ylim[2])
    
  }else{
    pMain <- ggplot(data = data, aes(x = data[ ,1], y=data[ ,2], colour = batch)) + 
      geom_point(shape = 16) + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
      ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
      scale_color_manual(values = color.mixo(1:10)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title)
    
    pTop <- ggplot(data, aes(x = data[ ,1], fill = batch)) + 
      geom_density(size = density.lwd, alpha=0.5) + ylab('Density') + 
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
            plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
            axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(xlim[1], xlim[2]) + labs(title = title)
    
    pRight <- ggplot(data, aes(x=data[ ,2], fill = batch)) + 
      geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() + ylab('Density') +
      theme(axis.title.x = element_text(size = rel(0.8)), 
            axis.title.y = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(ylim[1], ylim[2])
  }
  
  
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
                                legend.direction = 'vertical', 
                                legend.key.height = unit(0.2, 'cm'),
                                legend.key.width = unit(0.1, 'cm'),
                                legend.title = element_text(size = rel(legend.title.cex)),
                                legend.spacing.x = unit(0.1, 'cm'),
                                legend.spacing.y = unit(0.1, 'cm'),
                                legend.text = element_text(size = rel(legend.cex))))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  grid.arrange(pTop + theme(legend.position = 'none'), legend, pMain + 
                 theme(legend.position = 'none'), pRight + theme(legend.position = 'none'), 
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  
}



box_plot_fun = function(data = data, x = x, y = y, title = NULL, batch.legend.title = 'Batch', 
                        x.angle = 0, x.hjust = 0.5){
  ggplot(data = data, aes(x = x, y = y, fill = x)) + stat_boxplot(geom = "errorbar", width = 0.4) + 
    geom_boxplot() + scale_fill_manual(values = color.mixo(1:10)) + theme_bw() + 
    theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust), panel.grid = element_blank(),
          axis.title.x = element_blank(), axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,size = rel(1))) + 
    labs(fill = batch.legend.title, y = 'value',title = title) 
} 


RleMicroRna2 <- function (object, maintitle = NULL, batch = batch, xlab = NA,
                          legend = TRUE, cex.lab = 1.2, cex.xaxis = 1, 
                          cex.yaxis = 1, abline.lwd=0.5, legend.cex = 0.8,
                          xaxis.dist.ratio = 0.1, outcex = 1, title.cex = 1.3) 
{
  colorfill = color.mixo(batch)
  nARR = dim(object)[2]
  nGEN = dim(object)[1]
  y = apply(object, 1, median)
  mva = matrix(nrow = nGEN, ncol = nARR)
  for (i in 1:nARR) {
    x = object[, i]
    mva[ ,i] = (x - y)
  }
  med = apply(mva, 2, median)
  MIN = min(mva, na.rm = TRUE)
  MAX = max(mva, na.rm = TRUE)
  par(las = 3)
  plot(med, xlim = c(0, nARR + 1), ylim = c(MIN, MAX), axes = FALSE, 
       xlab = xlab, ylab = 'Deviations',cex.lab = cex.lab)
  colnames(mva) = colnames(object)
  res = boxplot(data.frame(mva), outline = TRUE, add = TRUE, col = colorfill,
                xaxt = 'n', outcex = outcex, cex.axis = cex.yaxis) #outcex for outlier
  axis(1, cex.axis = cex.xaxis, at = 1:ncol(object), labels = NA)
  points(med, type = 'p', col = 'blue', cex = outcex)
  lines(med, type = 'l', col = 'blue', lty = 'dotted')
  title(main = maintitle, cex.main = title.cex)
  abline(0, 0, col = 'red', lwd = abline.lwd)
  par(las = 0)
  end_point = 0.5 + ncol(object)  # add degrees to the x axis
  box.max = max(max(res$stats), max(res$out))
  box.min = min(min(res$stats), min(res$out))
  box.range = box.max - box.min
  text(seq(1.2, end_point, by = 1), par("usr")[3] - xaxis.dist.ratio*box.range, 
       srt = 60, adj = 1, xpd = TRUE,
       labels = paste(colnames(object)), cex = cex.xaxis)
  if(legend == TRUE){
    legend('topright', legend = unique(batch), pch=15, col = unique(colorfill), cex = legend.cex)
  }
}



percentileofscore = function(df, control.index){
  df.percentile = df
  df.percentile[1:nrow(df), 1:ncol(df)] = NA
  for(i in 1:ncol(df)){
    control = sort(df[control.index, i])
    for(j in 1:nrow(df)){
      percentile.strick = sum(control < df[j, i])/length(control)
      percentile.weak = (length(control) - sum(control > df[j, i]))/length(control)
      percentile = (percentile.strick + percentile.weak)/2
      df.percentile[j, i] = percentile
      
    }
  }
  return(df.percentile)
}

################

percentile_norm = function(data = data, batch = batch, trt = trt){
  batch = as.factor(batch)
  trt = as.factor(trt)
  trt.list = list()
  data.pn.df = data.frame()
  for(i in 1:nlevels(batch)){
    trt.each.b = trt[batch == levels(batch)[i]]
    trt.list[[i]] = trt.each.b
    data.each.b.pn = percentileofscore(data[batch == levels(batch)[i],], 
                                       which(trt.each.b == levels(trt.each.b)[1]))
    data.pn.df = rbind(data.pn.df,data.each.b.pn)
  }
  names(trt.list) = levels(batch)
  data.pn.df.reorder = data.pn.df[rownames(data), ]
  return(data.pn.df.reorder)
}

#Silhouette coefficient
# calculates silhouette width average
calc.sil = function(
    x, # the PC variates
    y1, y2 = NULL, # factor of interest, e.g. known batch info or known treatment info
    name.y1, name.y2 = NULL # character of the factor of interest
){
  library(cluster)
  # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
  dist.res = daisy(x, metric = 'euclidean')
  # for factor 1
  sil.batch.res1 = silhouette(x = as.numeric(y1), dist = dist.res)
  # if factor 2 is provided
  if(!is.null(y2)) sil.batch.res2 = silhouette(x = as.numeric(y2), dist = dist.res)
  
  # extract average width silhouette per level
  res1 = c(summary(sil.batch.res1)['clus.avg.widths']$clus.avg.widths)
  names(res1) = levels(y1)
  
  
  if(!is.null(y2)){
    res2 = c(summary(sil.batch.res2)['clus.avg.widths']$clus.avg.widths)
    names(res2) = levels(y2)
  }
  
  # output data for plotting
  if(!is.null(y2)){
    silh.coeff = c(res1, res2)
    Cluster = c(levels(y1), levels (y2))
    Type = c(rep(name.y1, nlevels(y1)), rep(name.y2, nlevels(y2)))
    data.plot = data.frame(silh.coeff, Cluster, Type)
    
  }else{
    silh.coeff = c(res1)
    Cluster = c(levels(y1))
    Type = rep(name.y1, nlevels(y1))
    data.plot = data.frame(silh.coeff, Cluster, Type)
  }
  
  return(invisible(data.plot))
}
### box&density function 
analyze_df <- function(df, plot_path, result_path) {
  
  color.mixo <- c('red', 'green', 'blue', 'yellow')
  
  # Check if the result path exists, create it if necessary
  if (!dir.exists(result_path)) {
    dir.create(result_path, recursive = TRUE)
  }
  
  for (col in 1:ncol(df)) {
    before.df <- data.frame(value = df[, col], batch = meta$EXP)
    box_plot_fun(data = before.df, x = before.df$batch,
                 y = before.df$value, title = paste('Column', col),
                 batch.legend.title = 'experiment')
    ggplot2::ggsave(paste(plot_path, '/box_plot_', col, '.png', sep = ''), width = 6, height = 4)
    
    ggplot(before.df, aes(x = value, fill = batch)) +
      geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo[1:10]) +
      labs(title = paste('Column', col), x = 'Value', fill = 'experiment') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                         panel.grid = element_blank())
    ggplot2::ggsave(paste(plot_path, '/density_plot_', col, '.png', sep = ''), width = 6, height = 4)
    
    data.lm <- lm(df[, col] ~ meta$EXP)
    summary_file <- paste(result_path, '/lm_summary_', col, '.txt', sep = '')
    sink(summary_file)
    summary(data.lm)
    sink(NULL)
  }
}


##INSTALL
BiocManager::install("")
##LIBS
library(knitr)
library(xtable) # table
library(mixOmics)
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(gridExtra) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(AgiMicroRna) # RLE plot
library(cluster) # silhouette coefficient
library(variancePartition) # variance calculation
library(pvca) # PVCA
library(pheatmap) # heatmap
library(ruv) # RUVIII
library(lmerTest) # lmer
library(bapred) # FAbatch
library(readxl)
library(MASS)
library(dabestr)
library(dplyr)
## read functions+libraries prior to analysis
##### real data
data = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/drp1_t.csv', header=TRUE)
meta = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/drp1_p.csv', header=TRUE)
##
data = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/mfn1_t.csv', header=TRUE)
meta = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/mfn1_p.csv', header=TRUE)
##
data = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/mfn2_t.csv', header=TRUE)
meta = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/mfn2_p.csv', header=TRUE)
##
data = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/vat_t.csv', header=TRUE)
meta = read.csv('D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/frames/vat_p.csv', header=TRUE)
## run 
type(data)
dim(data)
type(meta)
head(data)
head(meta)
rownames(data) <- data$sample
data = subset(data, select = -c(sample))

meta$EXP <- as.factor(meta$EXP)

data.clr <- logratio.transfo(data, logratio = 'CLR')
dim(data.clr)
class(data.clr) <- 'matrix' 
write.csv(as.data.frame(data.clr), file="data.clr-before-mfn2.csv")
#pca
data.pca.before <- pca(data.clr, ncomp = 2)

data.pca.plot.before <- Scatter_Density(data = data.pca.before$variates$X, batch = meta$EXP, 
                                        trt = meta$severity, expl.var = data.pca.before$explained_variance, 
                                        xlim = c(-4.5,5), ylim = c(-3,4), 
                                        batch.legend.title = 'EXP', 
                                        trt.legend.title = 'sample_type', 
                                        title = 'PCA_Before_batch_effect_correction (HD)')


###################################### box & density plots
#color load
color.mixo <- c('red', 'green', 'blue', 'yellow')
#define path and run
plot_path <- 'D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/plots_batch-MS-vat1'
result_path <- 'D:/MiguelW12/Documents/noam_may_analysis/noam_new_batch/res2_batch-exp-vat3/lm_summary_real'
analyze_df(data.clr, plot_path, result_path)

######################################
head(meta)
table(meta2$severity)
##RLE-
meta.trt <- meta2$severity
meta.trt <- as.factor(meta.trt)
meta.batch <- meta2$plate
meta.batch <- as.factor(meta.batch)

data.batch_hc <- meta.batch[meta.trt == 'HC']
data.batch_sev <- meta.batch[meta.trt == 'S'] 
data.batch_pre <- meta.batch[meta.trt == 'PRE'] 
data.batch_mild <- meta.batch[meta.trt == 'M'] 
#t if needed
data.batch_sevT <- meta.batch[meta.trt == 'ST'] 
data.batch_preT <- meta.batch[meta.trt == 'PRET'] 
data.batch_mildT <- meta.batch[meta.trt == 'MT'] 



data.before_hc <- data.clr[meta.trt == 'HC', ]
data.before_sev <- data.clr[meta.trt == 'S', ]
data.before_pre <- data.clr[meta.trt == 'PRE', ]
data.before_mild <- data.clr[meta.trt == 'M', ]
#t if needed
data.before_sevT <- data.clr[meta.trt == 'ST', ]
data.before_preT <- data.clr[meta.trt == 'PRET', ]
data.before_mildT <- data.clr[meta.trt == 'MT', ]



RleMicroRna2(object = t(data.before_hc), batch = data.batch_hc, 
             maintitle = 'noam_batch_analysis (type=HC)')


RleMicroRna2(object = t(data.before_sev), batch = data.batch_sev, 
             maintitle = 'noam_batch_analysis (type=sev)')

RleMicroRna2(object = t(data.before_pre), batch = data.batch_pre, 
             maintitle = 'noam_batch_analysis (type=pre)')

RleMicroRna2(object = t(data.before_mild), batch = data.batch_mild, 
             maintitle = 'noam_batch_analysis (type=mild)')
#t if needed
RleMicroRna2(object = t(data.before_sevT), batch = data.batch_sevT, 
             maintitle = 'noam_batch_analysis (type=sevT)')

RleMicroRna2(object = t(data.before_preT), batch = data.batch_preT, 
             maintitle = 'noam_batch_analysis (type=preT)')

RleMicroRna2(object = t(data.before_mildT), batch = data.batch_mildT, 
             maintitle = 'noam_batch_analysis (type=mildT)')
###########################################
table(meta$plate)
table(meta$severity)
##heatmap-  
# scale
maxs <- apply(data.clr, 2, max)    
mins <- apply(data.clr, 2, min)
data.clr.scale <- scale(data.clr, center = mins, scale = maxs - mins)
# scale on samples
maxs2 <- apply(t(data.clr.scale), 2, max)    
mins2 <- apply(t(data.clr.scale), 2, min)
data.clr.scale <- scale(t(data.clr.scale), center = mins2, scale = maxs2 - mins2)
# ano+colors
data.anno_col <- data.frame(Batch = meta.batch, Tissue = meta.trt)
data.anno_metabo_colors <- list(Batch = c( 'p1' = '#8b33f6', 'p2' = 'red',
                                          'p3' = 'deeppink3', 'p4' = "darkgreen", 'p5' = 'firebrick',
                                          'p6' = 'chartreuse', 'p7' = 'yellow',
                                          'p8' = '#00FF00', 'p9' = '#00B100', 'p10' = 'blue4'),
                                Tissue = c(HC = '#F0E442', S = '#33f6ed', PRE = 'darkseagreen', 
                                           M = 'darkorchid1', ST = "black", PRET = "#DF536B", MT = "gray62"))


#plot
pheatmap(data.clr.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 8,
         annotation_col = data.anno_col,
         annotation_colors = data.anno_metabo_colors,
         border_color = 'NA',
         main = 'HD data - Scaled')



########################################
# var calc 
head(data.info)
data.form <- ~ meta.trt + meta.batch
data.info <- as.data.frame(cbind(rownames(data.clr), meta.trt, meta.batch))
rownames(data.info) <- rownames(data.clr)

# before
data.varPart.before <- fitExtractVarPartModel(exprObj = t(data.clr), 
                                              formula = data.form, 
                                              data = data.info)
data.varmat.before <- as.matrix(data.varPart.before[ ,1:2])

data.variance <- c(as.vector(data.varmat.before))
data.variance <- cbind(variance = data.variance, 
                       Type = rep(c('group', 'EXP'), each = ncol(data.clr)),
                       method = rep(c('before'), each = 2*ncol(data.clr)))
# reorder levels  
data.variance <- as.data.frame(data.variance)
data.variance$method <- factor(data.variance$method, 
                                 levels = unique(data.variance$method))
data.variance$variance <- as.numeric(as.character(data.variance$variance))



ggplot(data.variance, aes(x = Type, y = variance, fill = Type)) + 
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.text = element_text(size = 12), panel.grid = element_blank(), 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) + 
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)

################################################# 
##########
# bci bootstrap plots 
#run per marker 
# UT markers
bci_df = read.csv('C:/Users/Uri8s/OneDrive/Documents/stats/noam_markers_stats/vat_ut_final.csv', header=TRUE)

two.group.unpaired <- 
  bci_df %>%
  dabest(index, drp1_intensity,  
         # The idx below passes "Control" as the control group, 
         # and "Group1" as the test group. The mean difference
         # will be computed as mean(Group1) - mean(Control1).
         idx = c("HCDMSO", "PREDMSO", "MDMSO", "SDMSO"), 
         paired = FALSE)

two.group.unpaired
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
two.group.unpaired.meandiff
## PLOT
plot(two.group.unpaired.meandiff, 
     color.column = index,
     rawplot.markersize = 3,
     rawplot.groupwidth = 0.6,
     rawplot.ylabel = "drp1 Intensity",
     effsize.ylabel = "Mean Difference",
     axes.title.fontsize = 16,
     palette = c("forestgreen","firebrick3", "dodgerblue3", "purple"),
     theme = ggplot2::theme_classic()
     
)

## T
bci_df = read.csv('C:/Users/Uri8s/OneDrive/Documents/stats/noam_markers_stats/drp1_T_final.csv', header=TRUE)

two.group.unpaired <- 
  bci_df %>%
  dabest(index, drp1_intensity,  
         # The idx below passes "Control" as the control group, 
         # and "Group1" as the test group. The mean difference
         # will be computed as mean(Group1) - mean(Control1).
         idx = c("HCDMSO", "PREDMSO", "PREmdivi1 25uM",
                 "MDMSO", "Mmdivi1 25uM",
                 "SDMSO", "Smdivi1 25uM"), 
         paired = FALSE)

two.group.unpaired
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
two.group.unpaired.meandiff

plot(two.group.unpaired.meandiff, 
     color.column = index,
     rawplot.markersize = 2,
     rawplot.groupwidth = 0.4,
     rawplot.ylabel = "drp1 Intensity",
     effsize.ylabel = "Mean Difference",
     tick.fontsize = 4,
     axes.title.fontsize = 16,
     palette = c("forestgreen", "purple", "dodgerblue3",
                 "firebrick3", "orchid",
                 
                 "steelblue1","indianred"),
     theme = ggplot2::theme_classic()
     
)
##########################################################################################





