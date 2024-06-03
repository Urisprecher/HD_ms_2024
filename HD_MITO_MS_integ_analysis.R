### HD integration-FIG4
gc()
##install
BiocManager::install("", force=TRUE)
##libraries
library(corrplot)
library(Hmisc)
library(psych)
library(factoextra)
library(AICcmodavg)
library(networkD3)
library(webshot)
library(webshot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(patchwork)
library(circlize)
library(DAAG)
library(coefplot)
library(performance)
library(visreg)
library(parameters)
library(gt)
library(see)
library(InformationValue)
library(randomForest)
library(rms)
library(nnet)
library(pheatmap)
library(creditmodel)
library(ggplot2)
library(networkD3)
library(webshot)
library(webshot2)
##read data
id_data_cor = read.csv('D:/MiguelW12/Documents/stats_script_try/HD_integ/final_data_integ_id_nohc.csv', header=TRUE)
head(id_data)
id_data = subset(id_data, select = -c(X))
head(id_data)
#cor
#extract numeric columns for correlation analysis
numeric_data <- id_data_cor[sapply(id_data_cor, is.numeric)]

## remove unwanted columns
numeric_data <- numeric_data[, apply(numeric_data, 2, sd) > 0]

###compute correlation 
correlation_matrix <- cor(numeric_data, method = "pearson")

####compute pearson p-values 
correlation_results <- corr.test(numeric_data, method = "pearson")
# save correlation matrix
write.csv(correlation_matrix, file = 'cor_mat_id-final.csv', row.names = TRUE)
# save p vals matrix
write.csv(correlation_results$p, file = 'cor_res_id-final.csv', row.names = TRUE)
##plot 
corrplot(correlation_matrix, method = "circle", type = "full", bg = "white", title = 'HD-cor_by_id',
         order = "hclust", hclust.method = 'single', addrect = 3, rect.col = 'black', rect.lwd = 5, 
         tl.col = "black", tl.cex = 1, tl.srt = 45, addCoef.col = "black", number.cex = 0.8, 
         p.mat = correlation_results$p, sig.level = 0.05, insig = "label_sig", pch = '*', pch.col= 'yellow', pch.cex = 2,
         
         )

######################################################################
## clust
clust_dat_id = read.csv('D:/MiguelW12/Documents/stats_script_try/HD_integ/clust-data_integ_id_nohc.csv', header=TRUE)
#PREP DATA
rownames(clust_dat_id) <- clust_dat_id$group_with_id
clust_dat_id = subset(clust_dat_id, select = -c(group_with_id))
df <- as.data.frame(clust_dat_id)
scaled_data <- scale(df)
#v1-no k 
clust <- pheatmap(scaled_data, cluster_cols = T, cluster_rows = T, main = "Cluster_Heatmap-HD", scale = "column",  
                  fontsize = 8, clustering_distance_rows = "canberra", clustering_distance_cols = "canberra", 
                  clustering_method = "single", show_colnames = TRUE, 
                  border_color = "black", cellwidth = 30, cellheight = 30, treeheight_row = 40, treeheight_col = 60, angle_col = 45, 
                  
                  )
#v2- with k
clust_k <- pheatmap(scaled_data, cluster_cols = T, cluster_rows = T, main = "Cluster_Heatmap-HD",  
                  fontsize = 8, clustering_distance_rows = "canberra", clustering_distance_cols = "canberra", 
                  clustering_method = "single", show_colnames = TRUE, kmeans_k = 3, cutree_rows = 2,
                  border_color = "black", treeheight_row = 40, treeheight_col = 60, angle_col = 45, 
                  
)
#save clust groups
clusts <- clust_k$kmeans$cluster
clusts
write.csv(as.data.frame(clusts), file = "hd_final-clustered_data$kmeans-clusts.csv", row.names = FALSE)




## SANKEY 


#read data
sank_data = read.csv('D:/MiguelW12/Documents/stats_script_try/HD_integ/sank_try.csv', header=TRUE)
#prep data
nodes4 <- data.frame(name=c(as.character(sank_data$source), as.character(sank_data$target)) %>% unique())
sank_data$IDsource=match(sank_data$source, nodes4$name)-1 
sank_data$IDtarget=match(sank_data$target, nodes4$name)-1
#color code
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF",
"#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
#plot
s <- sankeyNetwork(Links = sank_data, Nodes = nodes4,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)
s
#save html+pdf
saveNetwork(s, "sankey_diagram.html", selfcontained = FALSE)
webshot("sankey_diagram.html", file = "sankey_diagram.pdf", zoom = 2)


####################################################################################