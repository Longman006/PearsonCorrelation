library(ggplot2)
library(reshape2)

norm_vec <- function(x) sqrt(sum(x^2))

my.reshape<- function(frame){
  r.names<-row.names(frame)
  c.names<-colnames(frame)
  Var1<-array(0,length(r.names)*length(c.names))
  Var2<-array(0,length(r.names)*length(c.names))
  value<-array(0,length(r.names)*length(c.names))
  for (i in 1:length(r.names)) {
    for (j in 1:length(c.names)) {
      index<-(i*length(r.names))-1+j
      Var1[index] <- r.names[i]
      Var2[index] <- c.names[j]
      value[index] <- round(frame[i,j],2)
    }
  }
  return(data.frame(Var1 = Var1, Var2 = Var2, value = value))
}
Create.P.Matrix <- function(matrix, data){
  for(i in 1:nrow(data)){
    matrix[data$V1[i], data$V2[i]]=matrix[data$V1[i], data$V2[i]]+data$V4[i]
    matrix[data$V1[i], data$V3[i]]=matrix[data$V1[i], data$V3[i]]+data$V4[i]
  }
  return(matrix)
}

Pearson.Spearman.Correlation<-function(P.matrix,Rank){
  pearson.spearman<-array(0,c(1,nrow(P.matrix)-1))
  for(i in 1:(nrow(P.matrix)-1)){
    if(Rank == FALSE){
      a<-P.matrix[i,]-mean(P.matrix[i,])
      b<-P.matrix[i+1,]-mean(P.matrix[i+1,])
    }
    else{
      a<-rank(P.matrix[i,])-mean(rank(P.matrix[i,]))
      b<-rank(P.matrix[i+1,])-mean(rank(P.matrix[i+1,]))
    }
      pearson.spearman[i]<-t(a) %*% (b)/(norm_vec(a)*norm_vec(b))
  }
  
  return(pearson.spearman)
}

Pearson.Spearman.Correlation.Matrix<-function(P.matrix,Rank){
  pearson.spearman<-array(0,c(nrow(P.matrix),nrow(P.matrix)))
  for(i in 1:(nrow(P.matrix))){
    for(j in 1:(nrow(P.matrix))){
    if(Rank == FALSE){
      a<-P.matrix[i,]-mean(P.matrix[i,])
      b<-P.matrix[j,]-mean(P.matrix[j,])
    }
    else{
      a<-rank(P.matrix[i,])-mean(rank(P.matrix[i,]))
      b<-rank(P.matrix[j,])-mean(rank(P.matrix[j,]))
    }
    pearson.spearman[i,j]<-t(a) %*% (b)/(norm_vec(a)*norm_vec(b))
    }
  }
  
  return(pearson.spearman)
}

my.reorder <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return(cormat)
}
node.numbers<-c(6980,14488,16)
layer.sizes<-c(7,13,2)
layer.names<-list(arabidopsis_genetic_layers,arxiv_netscience_layers,Padgett.Florentine.Families_layers)
names<-c("Arabidopsis","Arxiv_netscience","Padgett")
data.all<-list(arabidopsis_genetic_multiplex,arxiv_netscience_multiplex,Padgett.Florentine.Families_multiplex)

P.matrices<-sapply(1:length(node.numbers), function(i) Create.P.Matrix(array(0,c(layer.sizes[i],node.numbers[i])),data.all[[i]]))

pearson.coefficients<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation(P.matrices[[i]],FALSE))
spearman.coefficients<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation(P.matrices[[i]],TRUE))

pearson.coefficients.matrix<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation.Matrix(P.matrices[[i]],FALSE))
spearman.coefficients.matrix<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation.Matrix(P.matrices[[i]],TRUE))

pearson.matrix<-sapply(1:length(pearson.coefficients.matrix),function(i) {
    data<-as.data.frame(pearson.coefficients.matrix[[i]])
    colnames(data)<-layer.names[[i]]$layerLabel
    rownames(data)<-layer.names[[i]]$layerLabel
    return(data)
})

melted.pearson.matrix <- lapply(1:length(pearson.matrix), function(i) my.reshape(pearson.matrix[[i]]))
ggheatmap <- ggplot(melted.pearson.matrix[[2]], aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    #legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

coefficients.list<-list(pearson.coefficients, spearman.coefficients)
coefficients.means<- sapply(coefficients.list, function(sp) sapply(sp, function(x) mean(x)))
df.means<-data.frame(name = c(replicate(3,"pearson"),replicate(3,"spearman")), coefficient = c(coefficients.means[,1], coefficients.means[,2]))


dfs<-lapply(1:3,function(i) data.frame(layer = 1:length(pearson.coefficients[[i]]), pearson = t(pearson.coefficients[[i]]),spearman = t(spearman.coefficients[[i]]),name = names[i]))
df<-rbind(dfs[[1]],dfs[[2]],dfs[[3]])

plots<-lapply(1:3,function(i){
  g<-ggplot(dfs[[i]], aes(x =layer )) +
  geom_point(aes(y = dfs[[i]]$pearson, color = "pearson"))+
  geom_line(aes(y = dfs[[i]]$pearson, color = "pearson"),alpha=0.25) +
  scale_colour_manual("", breaks = c("pearson", "spearman"), values = c("blue", "red"))+
  geom_smooth(aes(y = dfs[[i]]$pearson, color = "pearson"), fill="blue", alpha=0.2) +
  geom_point(aes(y = dfs[[i]]$spearman, color = "spearman"))+
  geom_line(aes(y = dfs[[i]]$spearman, color = "spearman"),alpha=0.25)+
  geom_smooth(aes(y = dfs[[i]]$spearman, color = "spearman"), fill="red", alpha=0.2)+
  labs(x = "Layer",y = "Coefficient")+
  ggtitle(dfs[[i]]$name)
  return(g)
})

boxplot.pearson<-ggplot(df, aes(x = name, y = pearson, color=name)) +
  geom_boxplot(outlier.shape=8, outlier.size=4)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=5)+
  geom_jitter(shape=16, position=position_jitter(0.2))

boxplot.spearman<-ggplot(df, aes(x = name, y = spearman, color=name)) +
  geom_boxplot(outlier.shape=8, outlier.size=4)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=5)+
  geom_jitter(shape=16, position=position_jitter(0.2))

boxplot.sum<-ggplot(df.means, aes(x = name, y = coefficient, color=name)) +
  geom_boxplot(outlier.shape=8, outlier.size=4)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=5)

all.boxplots<-list(boxplot.pearson,boxplot.spearman,boxplot.sum)



