library(cluster)
library(factoextra) 
library(dplyr)
library(ggplot2)
library(reshape2)

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

########################################################
# Cluster Analysis
########################################################

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

pdf(file = paste0(base, "average_silhouette coefficient.pdf"), width = 10, height = 7)
fviz_nbclust(df_meta_norm, kmeans, method = "silhouette")
dev.off()


k2 <- kmeans(df_meta_norm, centers = 5, nstart = 25)
str(k2)
pdf(file = paste0(base, "clusters.pdf"), width = 10, height = 10)
fviz_cluster(k2, data = df_meta_norm, geom = "point", ggtheme = theme_bw())
dev.off()

df_match[,'cluster'] <- k2$cluster
write.csv(df_match, file = paste0(base, "spot_matching_result_imputated_cluster.csv"))

##############################################################
#Evaluate clusters
##############################################################

# Number of spots for each cluster
cluster_histogram <- k2$cluster %>% table()
cluster_histogram_df <- data.frame(cluster = names(cluster_histogram), freq = as.numeric(cluster_histogram))

pdf(file = paste(base, "Cluster_freq_histogram.pdf", sep=""), width = 9, height = 7)
ggplot_cluster_histogram <- ggplot(cluster_histogram_df, aes(x=cluster, y=freq)) +
  geom_bar(stat="identity", width=1, color = "white") +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
  geom_text(aes(y = (freq+10), label = freq, fontface = 2), size=8) +
  labs(title = "Number of Spots in Each Cluster")

print(ggplot_cluster_histogram)
dev.off()

# Is there a certain pattern in cluster?
# Pie plot : the proportion of cluster in each channel
for (ch in c('C1', 'C2', 'C3', 'C4')){
  pdf(file = paste(base, paste0("cluster_proportion_in_channel_", ch, ".pdf"), sep=""), width = 9, height = 7)
  
  cluster <- df_match[!is.na(df_match[,paste0(ch, "_id")]),"cluster"]
  cluster_freq <- cluster %>% table()
  cluster_freq_df <- data.frame(cluster = names(cluster_freq), freq = as.numeric(cluster_freq))
  empty_cluster <- setdiff(c(1,2,3,4,5), names(cluster_freq))
  if (length(empty_cluster) > 0){
    cluster_freq_df[(length(names(cluster_freq))+1),"cluster"] <- empty_cluster
    cluster_freq_df[(length(names(cluster_freq))+1),"freq"] <- 0
    cluster_freq_df$cluster <- factor(cluster_freq_df$cluster, levels = c(1,2,3,4,5))
  }
  
  ggplot_cluster_proportion <- ggplot(cluster_freq_df, aes(x=cluster, y=freq)) +
    geom_bar(stat="identity", width=1, color = "white") +
    #scale_fill_manual(values = c("1" = "#F8766D", "2" = "#D39200", "3" = "#00BA38", "4" = "#00aedb", "5" = "#C77CFF")) + #"5" = "#00B9E3", , "7" = "#FF61C3"
    theme_bw() +
    theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
    geom_text(aes(y = (freq+10), label = freq, fontface = 2), size=8) +
    labs(title = paste0("Cluster Propotion in channel ", ch))
  
  print(ggplot_cluster_proportion)
  dev.off()
}
# Pie plot : the proportion of cluster in each n_match
pdf(file = paste(base, "cluster_proportion_in_n_match.pdf", sep=""), width = 9, height = 7)

n_match <- df_match[which(df_match$n_match > 1),c("cluster", "n_match")]
n_match_freq_df <- n_match %>% table() %>% melt() %>% data.frame()

n_match_freq_df$cluster <- factor(n_match_freq_df$cluster)
n_match_freq_df$n_match <- factor(n_match_freq_df$n_match)

ggplot_n_match_proportion <- ggplot(n_match_freq_df, aes(x=cluster, y=value, fill=cluster)) +
  geom_bar(position = "dodge", stat="identity", width=1, color = "white") +
  facet_wrap(~n_match) +
  scale_fill_manual(values = c("1" = "#F8766D", "2" = "#D39200", "3" = "#00BA38", "4" = "#00aedb", "5" = "#C77CFF")) + #"5" = "#00B9E3", , "7" = "#FF61C3"
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
  geom_text(aes(y = (value + 7), label = value, fontface = 2), size=3) +
  labs(title = "Cluster Propotion in n_match")

print(ggplot_n_match_proportion)
dev.off()

# Stacked bar plot : the proportion of cluster# in each channel combination 
# [C1, C2, C3, C4, (C1, C2), (C1, C3), (C1, C4), (C2, C3), (C2, C4), (C3, C4), (C1, C2, C3), (C1, C2, C4), (C1, C3, C4), (C2, C3, C4), (C1, C2, C3, C4)]


# Mean sigma & intensity values in each cluster

mean_df <- data.frame(C1_sigma = c(), C2_sigma = c(), C3_sigma = c(), C4_sigma = c(), C1_int = c(), C2_int = c(), C3_int = c(), C4_int = c())
for (c in 1:5){
  sigma <- df_match[which(df_match$cluster == c),c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma")]
  int <- df_match[which(df_match$cluster == c), c("C1_int", "C2_int", "C3_int", "C4_int")]
  
  mean_sigma <- apply(sigma, 2, mean)
  mean_int <- apply(int, 2, mean)
  
  for (d in 1:4){
    mean_df[c, names(mean_sigma)[d]] <- mean_sigma[d]
    mean_df[c, names(mean_int)[d]] <- mean_int[d]
  }
}
mean_df <- mean_df[,c(1,3,5,7,2,4,6,8)]
mean_df
