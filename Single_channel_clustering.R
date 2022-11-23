library(cluster)
library(factoextra) 
library(dplyr)
library(ggplot2)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
library(plyr)

#base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
#base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
base = "C:/Users/pc/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
vis_path = paste0(base, "analysis_20221121/")
file = "C2_result.csv"

#######################################################
# Explore data
#######################################################

df_single <- read.csv(paste0(base, file))
print(colnames(df_single))

# 1. Remove positions that the spot is out of range in all channels
df_single<- df_single[!(df_single$sigma..nm. < 20 | df_single$sigma..nm. > 250),]
df_single<- df_single[!(df_single$intensity..photon. <  500 | df_single$intensity..photon. > 100000),]

df_meta <- df_single[,c("intensity..photon.", "sigma..nm.")] #, "chi2", "uncertainty..nm.")]

# Pair plot
pdf(file = paste0(vis_path, "pairs_plot.pdf"), width = 7, height = 7)
pairs(df_meta)
dev.off()


# Compare Intensity/Sigma distribution in channels
sig_dist <- df_meta[,"sigma..nm."] %>% melt()

pdf(file = paste0(vis_path, "sigma_distribution_histogram.pdf"), width = 8, height = 6)
ggplot(sig_dist, aes(x = value)) +
  geom_histogram(binwidth=5, alpha = .5, position = "identity")+
  geom_vline(xintercept=mean(sig_dist$value),linetype = "dashed", size = 1) +
  xlab("sigma")
#geom_vline(data = CI_low, aes(xintercept=CI, colour = variable), size = 0.5) +
#geom_vline(data = CI_high, aes(xintercept=CI, colour = variable), size = 0.5)
dev.off()

pdf(file = paste0(vis_path, "sigma_distribution_barplot.pdf"), width = 8, height = 6)
ggplot(sig_dist, aes(y=value)) +
  geom_boxplot() +
  theme_bw()
dev.off()


int_dist <- df_meta[,"intensity..photon."] %>% melt()

pdf(file = paste0(vis_path, "intensity_distribution_histogram.pdf"), width = 8, height = 6)
ggplot(int_dist, aes(x = value)) +
  geom_histogram(binwidth=100, alpha = .5, position = "identity")+
  geom_vline(xintercept=mean(int_dist$value),linetype = "dashed", size = 1) +
  xlab("intensity")
dev.off()

pdf(file = paste0(vis_path, "intensity_distribution_barplot.pdf"), width = 8, height = 6)
ggplot(int_dist, aes(y=value)) +
  geom_boxplot() +
  theme_bw()
dev.off()

########################################################
# Cluster Analysis
########################################################

# Normalization
df_meta <- log10(df_meta)
means <- apply(df_meta,2,mean)
sds <- apply(df_meta,2,sd)
df_meta_norm <- scale(df_meta,center=means,scale=sds)

# Histogram: Sigma/intensity distribution after normalization
df_meta_norm_dist <- df_meta_norm %>% melt()
df_meta_norm_dist <- df_meta_norm_dist[,c(2,3)]
colnames(df_meta_norm_dist)[1] <- "variable"
#cdat <- ddply(df_meta_norm_dist, "variable", summarise, rating.mean = mean(value))

pdf(file = paste0(vis_path, "sigma_intensity_distribution_norm_histogram.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_dist, aes(x = value, fill = variable)) +
  geom_histogram(binwidth=.1, alpha = .5, position = "identity")
  #geom_vline(data = cdat, aes(xintercept=rating.mean, colour = variable),linetype = "dashed", size = 1) +
dev.off()

pdf(file = paste0(vis_path, "sigma_intensity_distribution_norm_barplot.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_dist, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()


# Calculate distance and determine cluster number
# Euclidean distance: observation with high values of features will be clustered together
# Correlation based distance: Identify clusters of observation regardless of their magnitude
distance = dist(df_meta_norm)
df_meta.hclust = hclust(distance, method = "ward.D2")

pdf(file = paste0(vis_path, "dendrogram.pdf"), width = 10, height = 7)
plot(df_meta.hclust,hang=-1, labels = FALSE, main='Hierarchical Cluster')
dev.off()

# Determine the best number of clusters
sil_coef_vector <- c()
for (num_cluster in 2:10){
  print(paste0("computing cluster number ", num_cluster))
  sil_cutree <- cutree(df_meta.hclust, k = num_cluster)
  sil_cl <- silhouette(sil_cutree ,distance)
  sil_cl_summary <- summary(sil_cl)
  sil_coefficient <- as.vector(sil_cl_summary$si.summary["Mean"])
  sil_coef_vector<- c(sil_coef_vector, sil_coefficient)
}

pdf(file = paste0(vis_path, "silhouette coefficient.pdf"), width = 10, height = 7)
plot(2:10,sil_coef_vector,"b", xlab="k", ylab="silhouette coefficient")
dev.off()

pdf(file = paste0(vis_path, "average_silhouette coefficient.pdf"), width = 10, height = 7)
fviz_nbclust(df_meta_norm, kmeans, method = "silhouette")
dev.off()

# K-means clustering
n_cluster = 2
k2 <- kmeans(df_meta_norm, centers = n_cluster, nstart = 25)  # centers = number of clusters
str(k2)

pdf(file = paste0(vis_path, "clusters.pdf"), width = 10, height = 10)
fviz_cluster(k2, data = df_meta_norm, geom = "point", ggtheme = theme_bw())
dev.off()

df_single[,'cluster'] <- factor(k2$cluster)

library(GGally)
pdf(file = paste0(vis_path, "clusters_pairs.pdf"), width = 10, height = 10)
ggpairs(df_single, columns=c("intensity..photon.", "sigma..nm."), aes(colour=cluster, alpha = 0.5), lower=list(continuous='points'), axisLabels='none')#, upper=list(continuous='blank'))
dev.off()

#write.csv(df_match, file = paste0(base, "spot_matching_result_imputated_cluster.csv"))

##############################################################
#Evaluate clusters
##############################################################

# Number of spots for each cluster
print(k2$size)

# Is there a certain pattern in cluster?
# Pie plot : the proportion of cluster in each n_match
cluster_histogram <- k2$cluster %>% table()
cluster_histogram_df <- data.frame(cluster = names(cluster_histogram), freq = as.numeric(cluster_histogram))

pdf(file = paste(vis_path, "cluster_proportion_in_n_match.pdf", sep=""), width = 9, height = 7)

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

# Mean sigma & intensity values in each cluster

mean_df <- data.frame(sigma = c(), int = c())
for (c in 1:n_cluster){
  sigma <- df_single[which(df_single$cluster == c),"sigma..nm."]
  int <- df_single[which(df_single$cluster == c),"intensity..photon."]
  
  mean_sigma <- mean(sigma)
  mean_int <- mean(int)
  
  mean_df[c, "sigma"] <- mean_sigma
  mean_df[c, "int"] <- mean_int
}
mean_df

