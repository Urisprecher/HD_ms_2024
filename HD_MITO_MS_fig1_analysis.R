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
#--------------------------
# function that calculates the silhouette coefficient based on a known cluster (i.e. batch or treatment)
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
## box&dist plot function 
analyze_df <- function(df, plot_path, result_path) {
  
  color.mixo <- c('red', 'green', 'blue', 'yellow')
  
  # Check if the result path exists, create it if necessary
  if (!dir.exists(result_path)) {
    dir.create(result_path, recursive = TRUE)
  }
  
  for (col in 1:ncol(df)) {
    before.df <- data.frame(value = df[, col], batch = meta$experiment)
    box_plot_fun(data = before.df, x = before.df$batch,
                 y = before.df$value, title = paste('Column', col),
                 batch.legend.title = 'experiment')
    ggsave(paste(plot_path, '/box_plot_', col, '.png', sep = ''), width = 6, height = 4)
    
    ggplot(before.df, aes(x = value, fill = batch)) +
      geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo[1:10]) +
      labs(title = paste('Column', col), x = 'Value', fill = 'experiment') +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                         panel.grid = element_blank())
    ggsave(paste(plot_path, '/density_plot_', col, '.png', sep = ''), width = 6, height = 4)
    
    data.lm <- lm(df[, col] ~ meta$experiment)
    summary_file <- paste(result_path, '/lm_summary_', col, '.txt', sep = '')
    sink(summary_file)
    summary(data.lm)
    sink(NULL)
  }
}

#install 
BiocManager::install("")
##required libs
library(Matrix)
library(knitr)
library(BayesFactor)
library(xtable) 
library(mixOmics)
library(sva) 
library(ggplot2) 
library(gridExtra) 
library(limma) 
library(vegan) 
library(AgiMicroRna) 
library(cluster) 
library(variancePartition)
library(pvca) 
library(pheatmap) 
library(ruv) 
library(lmerTest) 
library(bapred) 
library(readxl)
library(dabestr)
library(dplyr)
## read functions+libraries prior to analysis
##### read data
data = read.csv('D:/MiguelW12/Documents/Adam/adam_batch/adam_mdivi_t.csv', header=TRUE)
meta = read.csv('D:/MiguelW12/Documents/Adam/adam_batch/adam_mdivi_p.csv', header=TRUE)
type(data)
dim(data)
type(meta)
dim(meta)

rownames(data) <- data$sample
data = subset(data, select = -c(sample))


data.clr <- logratio.transfo(data, logratio = 'CLR')
class(data.clr) <- 'matrix' 
write.csv(as.data.frame(data.clr), file="adam_mdivi-data.clr-before.csv")
meta$stage
meta$experiment
data.pca.before <- pca(data.clr, ncomp = 3)

data.pca.plot.before <- Scatter_Density(data = data.pca.before$variates$X, batch = meta$experiment, 
                                        trt = meta$stage, expl.var = data.pca.before$explained_variance, 
                                        xlim = c(-4.5,5), ylim = c(-3,4), 
                                        batch.legend.title = 'EXP', 
                                        trt.legend.title = 'sample_type', 
                                        title = 'PCA_Before_batch_effect_correction (HD)')


######################################
## box plots + dist plots for each feature
#color selection
color.mixo <- c('red', 'green', 'blue', 'yellow')
#define path and run
plot_path <- 'D:/MiguelW12/Documents/Adam/adam_may_analysis/batch/plots'
result_path <- 'D:/MiguelW12/Documents/Adam/adam_may_analysis/batch/res2/lm_summary'
analyze_df(data.clr, plot_path, result_path)
######################################
##RLE- 
meta.trt <- meta$stage
meta.trt <- as.factor(meta.trt)
meta.batch <- meta$experiment
meta.batch <- as.factor(meta.batch)
table(meta.trt)
table(meta.batch)

data.batch_hc <- meta.batch[meta.trt == 'HC']
data.batch_sev <- meta.batch[meta.trt == 'Severe'] 
data.batch_pre <- meta.batch[meta.trt == 'Premanifest'] 
data.batch_mild <- meta.batch[meta.trt == 'Mild'] 



data.before_hc <- data.clr[meta.trt == 'HC', ]
data.before_sev <- data.clr[meta.trt == 'Severe', ]
data.before_pre <- data.clr[meta.trt == 'Premanifest', ]
data.before_mild <- data.clr[meta.trt == 'Mild', ]


RleMicroRna2(object = t(data.before_hc), batch = data.batch_hc, 
             maintitle = 'adam_batch_analysis (type=HC)')

RleMicroRna2(object = t(data.before_sev), batch = data.batch_sev, 
             maintitle = 'adam_batch_analysis (type=sev)')

RleMicroRna2(object = t(data.before_pre), batch = data.batch_pre, 
             maintitle = 'adam_batch_analysis (type=pre)')

RleMicroRna2(object = t(data.before_mild), batch = data.batch_mild, 
             maintitle = 'adam_batch_analysis (type=mild)')


###########################################


######################################## var calculations

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
                       Type = rep(c('sample', 'experiment'), each = ncol(data.clr)),
                       method = rep(c('Before'), each = 2*ncol(data.clr)))
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

#######################################
######## BOOTSTRAP PLOTS 
# PER FEATURE 
#ADDbci
bci_df = read.csv('C:/Users/Uri8s/OneDrive/Documents/stats/noam_markers_stats/adam_mdivi_stats.csv', header=TRUE)
#### run per feature
# bootstrap plots 
two.group.unpaired <- 
  bci_df %>%
  dabest(index, network_count_per_cell,  
         # The idx below passes "Control" as the control group, 
         # and "Group1" as the test group. The mean difference
         # will be computed as mean(Group1) - mean(Control1).
         idx = c("HCDMSO", 
                 "PremanifestDMSO", 
                 "MildDMSO", 
                 "SevereDMSO" ), 
         paired = FALSE)

two.group.unpaired
two.group.unpaired.meandiff <- mean_diff(two.group.unpaired)
two.group.unpaired.meandiff

## final plot

plot(two.group.unpaired.meandiff, 
     color.column = index,
     rawplot.markersize = 3,
     rawplot.groupwidth = 0.6,
     rawplot.ylabel = "network_count_per_cell",
     effsize.ylabel = "Mean Difference",
     axes.title.fontsize = 16,
     palette = c("dodgerblue3","forestgreen", "purple", "firebrick3"),
     theme = ggplot2::theme_classic()
     
)


####################################################

