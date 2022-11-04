library(cluster)
library(factoextra) 

base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
file = "Spot_matching_result_imputated.csv"

df_match <- read.csv(paste0(base, file))
print(colnames(df_match))

df_meta <- df_match[,c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma", "C1_int", "C2_int", "C3_int", "C4_int")]

pdf(file = paste0(base, "pairs_plot.pdf"), width = 7, height = 7)
pairs(df_meta)
dev.off()

df_meta <- scale(df_meta)
head(df_meta)

# Normalization
means <- apply(df_meta,2,mean)
sds <- apply(df_meta,2,sd)
df_meta_norm <- scale(df_meta,center=means,scale=sds)

# Calculate distance and determine cluster number
distance = dist(df_meta_norm)
df_meta.hclust = hclust(distance)

pdf(file = paste0(base, "dendrogram.pdf"), width = 10, height = 7)
plot(df_meta.hclust,hang=-1, labels = FALSE, main='Hierarchical Cluster')
dev.off()

# Determine the best number of clusters
sil_coef_vector <- c()
for (num_cluster in 2:10){
  sil_cutree <- cutree(df_meta.hclust, k = num_cluster)
  sil_cl <- silhouette(sil_cutree ,distance)
  sil_cl_summary <- summary(sil_cl)
  sil_coefficient <- as.vector(sil_cl_summary$si.summary["Mean"])
  sil_coef_vector<- c(sil_coef_vector, sil_coefficient)
}

pdf(file = paste0(base, "silhouette coefficient.pdf"), width = 10, height = 7)
plot(2:10,sil_coef_vector,"b", xlab="k", ylab="silhouette coefficient")
dev.off()


fviz_nbclust(df_meta_norm, kmeans, method = "silhouette")

# Cluster analysis

k2 <- kmeans(df_meta_norm, centers = 5, nstart = 25)
str(k2)
pdf(file = paste0(base, "clusters.pdf"), width = 10, height = 10)
fviz_cluster(k2, data = df_meta_norm, geom = "point", ggtheme = theme_bw())
dev.off()

df_match[,'cluster'] <- k2$cluster
write.csv(df_match, file = paste0(base, "spot_matching_result_imputated_cluster.csv"))

# Evaluate clusters

# Is there a certain pattern in cluster?
# Pie plot : the proportion of cluster# in each channel
# Pie plot : the proportion of n_match in each cluster#
# Stacked bar plot : the proportion of cluster# in each channel combination 
# [C1, C2, C3, C4, (C1, C2), (C1, C3), (C1, C4), (C2, C3), (C2, C4), (C3, C4), (C1, C2, C3), (C1, C2, C4), (C1, C3, C4), (C2, C3, C4), (C1, C2, C3, C4)]
# Mean sigma & intensity values in each cluster