
library(Mfuzz)


metabolite <- read.csv("", row.names = 1, check.names = FALSE)
metabolite <- as.matrix(metabolite)



## Filtering----
#
mfuzz_class <- new('ExpressionSet',exprs = metabolite)

mfuzz_class <- filter.std(mfuzz_class,min.std=0)




mfuzz_class <- standardise(mfuzz_class)


set.seed(123)
cluster_num <- 4
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))




mfuzz.plot2(mfuzz_class,lwd=2, xlab="Stage", mfrow=c(2,2), ylab="Chemical changes",cl = mfuzz_cluster,time.labels = colnames(metabolite),
            centre=T, centre.col="black",centre.lwd=2)



cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size


head(mfuzz_cluster$cluster)


head(mfuzz_cluster$membership)


protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(metabolite[names(protein_cluster), ], protein_cluster)
head(protein_cluster)
write.csv(protein_cluster, '')




protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
head(protein_standard_cluster)
write.csv(protein_standard_cluster, 'evolution_standarized_cluster.csv')
##--------------------------------------------------------------------------------------------------------------------

## plot standarized cluster 1 and cluster4
cluster1 <- read.csv("",check.names = FALSE)
cluster1<-cluster1[,-c(2,8)]
cluster2 <- read.csv("",check.names = FALSE)
cluster2<-cluster2[,-c(2,3,9)]
cluster3 <- read.csv("",check.names = FALSE)
cluster3<-cluster3[,-c(2,3,9)]
cluster4 <- read.csv("",check.names = FALSE)
cluster4<-cluster4[,-c(2,8)]

df_cluster1=reshape2::melt(cluster1,
                           id.vars=1, 
                  measure.vars=2:6)
df_cluster2=reshape2::melt(cluster2,
                           id.vars=1, 
                           measure.vars=2:6)
df_cluster3=reshape2::melt(cluster3,
                           id.vars=1, 
                           measure.vars=2:6)
df_cluster4=reshape2::melt(cluster4,
                           id.vars=1, 
                           measure.vars=2:6)

library(ggplot2)
p1=ggplot(data=df_cluster1,
         aes(x=variable,y=value,
             group=chemicals,colour=chemicals))+
  geom_point(size=3)+
  labs(x="stage", y="z score")+
  geom_line()+
  scale_x_discrete(limits=c("HC","EPL","EE","AE","IE"))+
  theme_bw() 
p1


p2=ggplot(data=df_cluster2,
          aes(x=variable,y=value,
              group=chemicals,colour=chemicals))+
  geom_point(size=3)+
  labs(x="stage", y="z score")+
  geom_line()+
  scale_x_discrete(limits=c("HC","EPL","EE","AE","IE"))+
  theme_bw()
p2



p3=ggplot(data=df_cluster3,
          aes(x=variable,y=value,
              group=chemicals,colour=chemicals))+
  geom_point(size=3)+
  labs(x="stage", y="z score")+
  geom_line()+
  scale_x_discrete(limits=c("HC","EPL","EE","AE","IE"))+
  theme_bw()
p3


p4=ggplot(data=df_cluster4,
          aes(x=variable,y=value,
              group=chemicals,colour=chemicals))+
  geom_point(size=3)+
  labs(x="stage", y="z score")+
  geom_line()+
  scale_x_discrete(limits=c("HC","EPL","EE","AE","IE"))+
  theme_bw()
p4


