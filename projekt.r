library(ggplot2)

norm_vec <- function(x) sqrt(sum(x^2))

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

node.numbers<-c(6980,14488,16)
layer.sizes<-c(7,13,2)
data.all<-list(arabidopsis_genetic_multiplex,arxiv_netscience_multiplex,Padgett.Florentine.Families_multiplex)

P.matrices<-sapply(1:length(node.numbers), function(i) Create.P.Matrix(array(0,c(layer.sizes[i],node.numbers[i])),data.all[[i]]))

pearson.coefficients<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation(P.matrices[[i]],FALSE))
spearman.coefficients<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation(P.matrices[[i]],TRUE))

pearson.coefficients.matrix<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation.Matrix(P.matrices[[i]],FALSE))
spearman.coefficients.matrix<-sapply(1:length(P.matrices),function(i) Pearson.Spearman.Correlation.Matrix(P.matrices[[i]],TRUE))

coefficients.list<-list(pearson.coefficients, spearman.coefficients)
coefficients.means<- sapply(coefficients.list, function(sp) sapply(sp, function(x) mean(x)))
df.means<-data.frame(name = c(replicate(3,"pearson"),replicate(3,"spearman")), coefficient = c(coefficients.means[,1], coefficients.means[,2]))

names<-c("Arabidopsis","Arxiv_netscience","Padgett")
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



