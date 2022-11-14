library(cluster)
library(factoextra) 
library(dplyr)
library(ggplot2)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)
library(plyr)

#base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Negative/2/"
base = "C:/Users/user/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
# = "C:/Users/pc/Desktop/UNIST_internship/Sample_Image/Positive/PB417_01/"
vis_path = paste0(base, "analysis_20221114_MixtureModel/")
file = "Spot_matching_result_imputated.csv"

#######################################################
# Explore data
#######################################################

df_match <- read.csv(paste0(base, file))
print(colnames(df_match))

# Limit Ranges:
# 1. Limit sigma (30<sigma<250)
index_C1 <- which(df_match$C1_sigma < 30 | df_match$C1_sigma > 250)
df_match$C1_sigma[index_C1] <- 0
df_match$C1_int[index_C1] <- 0

index_C2 <- which(df_match$C2_sigma < 30 | df_match$C2_sigma > 250)
df_match$C2_sigma[index_C2] <- 0
df_match$C2_int[index_C2] <- 0

index_C3 <- which(df_match$C3_sigma < 30 | df_match$C3_sigma > 250)
df_match$C3_sigma[index_C3] <- 0
df_match$C3_int[index_C3] <- 0

index_C4 <- which(df_match$C4_sigma < 30 | df_match$C4_sigma > 250)
df_match$C4_sigma[index_C4] <- 0
df_match$C4_int[index_C4] <- 0

# 2. Remove positions that the spot is out of range in all channels
df_match <- df_match[!(df_match$C1_sigma==0 & df_match$C2_sigma==0 & df_match$C3_sigma==0 & df_match$C4_sigma==0),]

#df_meta <- df_match[,c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma", "C1_int", "C2_int", "C3_int", "C4_int")]
df_meta <- df_match[,c("C1_int", "C2_int", "C3_int", "C4_int")]
# Pair plot
pdf(file = paste0(vis_path, "pairs_plot.pdf"), width = 7, height = 7)
pairs(df_meta)
dev.off()

# Venn Diagram
C1 <- df_match[!is.na(df_match$C1_id), "X"]
C2 <- df_match[!is.na(df_match$C2_id), "X"]
C3 <- df_match[!is.na(df_match$C3_id), "X"]
C4 <- df_match[!is.na(df_match$C4_id), "X"]

myCol <- brewer.pal(4, "Set1")

# Chart
venn.diagram(
  x = list(C1, C2, C3, C4),
  category.names = c("C1" , "C2 " , "C3", "C4"),
  filename = paste0(vis_path, 'venn_diagram.png'),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
)


# Number of spots in each channel
n_spots <- data.frame(channel = c('C1', 'C2', 'C3', 'C4'), freq = c(length(C1), length(C2), length(C3), length(C4)))

pdf(file = paste0(vis_path, "n_spots_each_channel.pdf"), width = 8, height = 6)
ggplot(n_spots, aes(x = channel, y = freq)) + 
  geom_bar(position = "dodge", stat="identity", width=1, color = "white") +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
  geom_text(aes(y = (freq + 1000), label = freq, fontface = 2), size=5) +
  labs(title = "Number of Spots in Each Channel")
dev.off()


# Average intensity/sigma of spots in each channel
avg_int_prep <- df_match[,c("C1_int", "C2_int", "C3_int", "C4_int")]
avg_sig_prep <- df_match[,c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma")]
avg_sig_prep[avg_sig_prep == 0] <- NA

avg_int <- apply(avg_int_prep, 2, mean, na.rm = TRUE)
avg_sig <- apply(avg_sig_prep, 2, mean, na.rm = TRUE)

df_avg_int <- data.frame(Channel = c('C1', 'C2', 'C3', 'C4'), Average_Intensity = avg_int)
df_avg_sig <- data.frame(Channel = c('C1', 'C2', 'C3', 'C4'), Average_Sigma = avg_sig)

df_avg_int$Average_Intensity <- round(df_avg_int$Average_Intensity,2)
df_avg_sig$Average_Sigma <- round(df_avg_sig$Average_Sigma,2)

pdf(file = paste0(vis_path, "avg_int_histogram.pdf"), width = 8, height = 6)
ggplot(df_avg_int, aes(x = Channel, y = Average_Intensity)) + 
  geom_bar(position = "dodge", stat="identity", width=1, color = "white") +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
  geom_text(aes(y = (Average_Intensity + 1000), label = Average_Intensity, fontface = 2), size=5) +
  labs(title = "Average Intensity in Each Channel")
dev.off()

pdf(file = paste0(vis_path, "avg_sig_histogram.pdf"), width = 8, height = 6)
ggplot(df_avg_sig, aes(x = Channel, y = Average_Sigma)) + 
  geom_bar(position = "dodge", stat="identity", width=1, color = "white") +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
  geom_text(aes(y = (Average_Sigma + 3), label = Average_Sigma, fontface = 2), size=5) +
  labs(title = "Average Sigma in Each Channel")
dev.off()


# Compare Intensity/Sigma distribution in channels
multi_sig_dist <- avg_sig_prep %>% melt()
multi_sig_dist <- multi_sig_dist[which(multi_sig_dist$value > 0),]
cdat <- ddply(multi_sig_dist, "variable", summarise, rating.mean = mean(value))
#CI_low <- ddply(multi_sig_dist, "variable", summarise, CI = as.numeric(t.test(value, conf.level = 0.95)[4]$conf.int[1]))
#CI_high <- ddply(multi_sig_dist, "variable", summarise, CI = as.numeric(t.test(value, conf.level = 0.95)[4]$conf.int[2]))

pdf(file = paste0(vis_path, "sigma_distribution_histogram.pdf"), width = 8, height = 6)
ggplot(multi_sig_dist, aes(x = value, fill = variable)) +
  geom_histogram(binwidth=5, alpha = .5, position = "identity")+
  geom_vline(data = cdat, aes(xintercept=rating.mean, colour = variable),linetype = "dashed", size = 1) +
  xlab("sigma")
  #geom_vline(data = CI_low, aes(xintercept=CI, colour = variable), size = 0.5) +
  #geom_vline(data = CI_high, aes(xintercept=CI, colour = variable), size = 0.5)
dev.off()

pdf(file = paste0(vis_path, "sigma_distribution_barplot.pdf"), width = 8, height = 6)
ggplot(multi_sig_dist, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()

multi_int_dist <- avg_int_prep %>% melt()
multi_int_dist <- multi_int_dist[which(multi_int_dist$value > 0),]
cdat <- ddply(multi_int_dist, "variable", summarise, rating.mean = mean(value))

pdf(file = paste0(vis_path, "intensity_distribution_histogram.pdf"), width = 8, height = 6)
ggplot(multi_int_dist, aes(x = value, fill = variable)) +
  geom_histogram(binwidth=1000, alpha = .5, position = "identity")+
  geom_vline(data = cdat, aes(xintercept=rating.mean, colour = variable),linetype = "dashed", size = 1) +
  xlim(0, 200000) +
  ylim(0, 150) +
  xlab("intensity")
dev.off()

pdf(file = paste0(vis_path, "intensity_distribution_barplot.pdf"), width = 8, height = 6)
ggplot(multi_int_dist, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() +
  ylim(0, 200000) +
  theme_bw()
dev.off()

########################################################
# Cluster Analysis
########################################################

# Normalization
means <- apply(df_meta,2,mean)
sds <- apply(df_meta,2,sd)
df_meta_norm <- scale(df_meta,center=means,scale=sds)

# Histogram: Sigma/intensity distribution after normalization
df_meta_norm_sig <- df_meta_norm[,c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma")]
df_meta_norm_sig_dist <- df_meta_norm_sig %>% melt()
df_meta_norm_sig_dist <- df_meta_norm_sig_dist[,c(2,3)]
colnames(df_meta_norm_sig_dist)[1] <- "variable"
df_meta_norm_sig_dist <- df_meta_norm_sig_dist[c(rownames(multi_sig_dist)),]
cdat <- ddply(df_meta_norm_sig_dist, "variable", summarise, rating.mean = mean(value))

pdf(file = paste0(vis_path, "sigma_distribution_norm_histogram.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_sig_dist, aes(x = value, fill = variable)) +
  geom_histogram(binwidth=.1, alpha = .5, position = "identity")+
  geom_vline(data = cdat, aes(xintercept=rating.mean, colour = variable),linetype = "dashed", size = 1) +
  xlab("sigma") +
  xlim(0, 10)
dev.off()

pdf(file = paste0(vis_path, "sigma_distribution_norm_barplot.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_sig_dist, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()

df_meta_norm_int <- df_meta_norm[,c("C1_int", "C2_int", "C3_int", "C4_int")]
df_meta_norm_int_dist <- df_meta_norm_int %>% melt()
df_meta_norm_int_dist <- df_meta_norm_int_dist[,c(2,3)]
colnames(df_meta_norm_int_dist)[1] <- "variable"
df_meta_norm_int_dist <- df_meta_norm_int_dist[c(rownames(multi_int_dist)),]
cdat <- ddply(df_meta_norm_int_dist, "variable", summarise, rating.mean = mean(value))

pdf(file = paste0(vis_path, "intensity_distribution_norm_histogram.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_int_dist, aes(x = value, fill = variable)) +
  geom_histogram(binwidth=.1, alpha = .5, position = "identity")+
  geom_vline(data = cdat, aes(xintercept=rating.mean, colour = variable),linetype = "dashed", size = 1) +
  xlab("sigma") +
  xlim(0, 15)
dev.off()

pdf(file = paste0(vis_path, "intensity_distribution_norm_barplot.pdf"), width = 8, height = 6)
ggplot(df_meta_norm_int_dist, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() +
  theme_bw() +
  ylim(0, 20)
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

# # Fuzzy c-means clustering
# 
# n_cluster = 2
# k2.fcm <- fcm(df_meta_norm, centers=n_cluster, nstart=5)
# plotcluster(k2.fcm, cp=1, trans=TRUE)
# 
# k2 <- ppclust2(k2.fcm, "kmeans")
# fviz_cluster(k2, data = df_meta_norm, 
#                          ellipse.type = "convex",
#                          palette = "jco",
#                          geom = "point")

df_match[,'cluster'] <- factor(k2$cluster)

library(GGally)
pdf(file = paste0(vis_path, "clusters_pairs.pdf"), width = 10, height = 10)
#ggpairs(df_match, columns=c("C1_sigma", "C2_sigma", "C3_sigma", "C4_sigma", "C1_int", "C2_int", "C3_int", "C4_int"), aes(colour=cluster, alpha = 0.5), lower=list(continuous='points'), axisLabels='none')#, upper=list(continuous='blank'))
ggpairs(df_match, columns=c("C1_int", "C2_int", "C3_int", "C4_int"), aes(colour=cluster, alpha = 0.5), lower=list(continuous='points'), axisLabels='none')#, upper=list(continuous='blank'))

dev.off()

write.csv(df_match, file = paste0(base, "spot_matching_result_imputated_cluster.csv"))

##############################################################
#Evaluate clusters
##############################################################

# Number of spots for each cluster
cluster_histogram <- k2$cluster %>% table()
cluster_histogram_df <- data.frame(cluster = names(cluster_histogram), freq = as.numeric(cluster_histogram))

pdf(file = paste(vis_path, "Cluster_freq_histogram.pdf", sep=""), width = 9, height = 7)
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
  pdf(file = paste(vis_path, paste0("cluster_proportion_in_channel_", ch, ".pdf"), sep=""), width = 9, height = 7)
  
  cluster <- df_match[!is.na(df_match[,paste0(ch, "_id")]),"cluster"]
  cluster_freq <- cluster %>% table()
  cluster_freq_df <- data.frame(cluster = names(cluster_freq), freq = as.numeric(cluster_freq))
  empty_cluster <- setdiff(c(1:n_cluster), names(cluster_freq))
  if (length(empty_cluster) > 0){
    cluster_freq_df[(length(names(cluster_freq))+1),"cluster"] <- empty_cluster
    cluster_freq_df[(length(names(cluster_freq))+1),"freq"] <- 0
    cluster_freq_df$cluster <- factor(cluster_freq_df$cluster, levels = c(1:n_cluster))
  }
  
  ggplot_cluster_proportion <- ggplot(cluster_freq_df, aes(x=cluster, y=freq)) +
    geom_bar(stat="identity", width=1, color = "white") +
    #scale_fill_manual(values = c("1" = "#F8766D", "2" = "#D39200", "3" = "#00BA38", "4" = "#00aedb", "5" = "#C77CFF")) + #"5" = "#00B9E3", , "7" = "#FF61C3"
    theme_bw() +
    theme(plot.title = element_text(size=20, hjust=0.5, face = "bold")) +
    geom_text(aes(y = (freq+10), label = freq, fontface = 2), size=8) +
    labs(title = paste0("Cluster Proportion in channel ", ch))
  
  print(ggplot_cluster_proportion)
  dev.off()
}
# Pie plot : the proportion of cluster in each n_match
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

# Stacked bar plot : the proportion of cluster# in each channel combination 
# [C1, C2, C3, C4, (C1, C2), (C1, C3), (C1, C4), (C2, C3), (C2, C4), (C3, C4), (C1, C2, C3), (C1, C2, C4), (C1, C3, C4), (C2, C3, C4), (C1, C2, C3, C4)]


# Mean sigma & intensity values in each cluster

mean_df <- data.frame(C1_sigma = c(), C2_sigma = c(), C3_sigma = c(), C4_sigma = c(), C1_int = c(), C2_int = c(), C3_int = c(), C4_int = c())
for (c in 1:n_cluster){
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

