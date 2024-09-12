# for tissue


# GGM (Gaussian graphical modeling.)

library("GeneNet")

Tissue.data <- read.csv("", row.names = 1,check.names = FALSE)

Tissue.data<- as.matrix(Tissue.data)

# Compute GGM 
tissue.pcor <- ggm.estimate.pcor(Tissue.data,method="static")

#estimating optimal shrinkage
met.edges<-network.test.edges(tissue.pcor,direct=TRUE)

write.csv(met.edges,"tissue_correlation.csv")
#Extract network containing significant association (i.e edges) with probability>0.9 (i.e. FDR<0.1)
# net<-extract.network(met.edges, cutoff.ggm = 0.9) 


#--------------------------------------------------------------------------------------------------------
## plot for tissue

library(dplyr)
library(tidyr)
library(ggplot2)




diameter.1 <- 8
nodes.1 <- read.csv("") %>%  
filter(Part != "Inorganic") %>% 
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, 5 * pi / 2 - (2 * pi / nrow(.)), length.out = nrow(.))) %>%
  mutate(x = (diameter.1/2) * cos(angle), 
         y = (diameter.1/2) * sin(angle))



diameter.2 <- 2   
nodes.2 <- read.csv("") %>%
  filter(Part == "Inorganic") %>%
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, 5 * pi / 2 - (2 * pi / nrow(.)), length.out = nrow(.))) %>%
  mutate(x = (diameter.2/2) * cos(angle) + 7, 
         y = (diameter.2/2) * sin(angle) + 2.25)


data.node.plot <- rbind(nodes.1, nodes.2)[-1] 


edges <- read.csv("", row.names = 1) 

edges <- edges %>%
  mutate(From.x = data.node.plot$x[match(edges$source, data.node.plot$label)],
         From.y = data.node.plot$y[match(edges$source, data.node.plot$label)],
         To.x = data.node.plot$x[match(edges$target, data.node.plot$label)],
         To.y = data.node.plot$y[match(edges$target, data.node.plot$label)]) 


#connected_nodes <- unique(c(edges$source, edges$target))


#data.node.plot.connected <- data.node.plot[data.node.plot$label %in% connected_nodes, ]


ggplot(data = data.node.plot) +
  geom_segment(data = edges, 
               aes(xend= To.x, yend = To.y, x = From.x, y = From.y, color = factor(correlation), linewidth =factor(correlation)),
               show.legend = F)+
  scale_linewidth_discrete(range = c(0.07, 0.25))+
  geom_point(aes(x = x, y = y, fill = category),
             show.legend = F, size=4, shape = 21, stroke = 0)+ 
  coord_fixed()+
  scale_color_manual(values = c("Negative" = "#427AB2", "Positive" = "#DBDB8D")) + 
  scale_fill_manual(values = c("exogenous compound" = "#b5b5d8",
                               "endogenous compound" = "#f6f2c5", 
                               "NEMs" = "#ffcca6", 
                               "EMs" = "#dbcfe4"))+  
  #geom_text(data = data.frame(part = c("Organic", "Inorganic"),
                              #x = c(0,7),
                             # y = c(0,2.25)),
            #aes(label = part, x = x, y = y), size = 3)+ 
  ##geom_text(data = data.frame(label = c( paste("Nodes:", nrow(data.node.plot)),paste("Edges:", nrow(edges)), #统计节点、连线),
    ## x = 4,
    ##y = c(-3,-3.5,-4)),
   ## aes(x = x, y = y, label = label),hjust = 0, size = 4.1)+
  theme_void() 

#
#
#
#
ggsave("Tissue_PLOT2.pdf", width = 5.5, height = 3)

ggsave("Tissue_PLOT2.png", width = 5.5, height = 3)



##---------------------------------------------------------------------------------------------------------------

# for plasma

# GGM (Gaussian graphical modeling.)

library("GeneNet")

Plasma.data <- read.csv("", row.names = 1,check.names = FALSE)

Plasma.data<- as.matrix(Plasma.data)

# Compute GGM 
Plasma.pcor <- ggm.estimate.pcor(Plasma.data,method="static")

#estimating optimal shrinkage
met.edges.plasma<-network.test.edges(Plasma.pcor,direct=TRUE)

write.csv(met.edges.plasma,"Plasma_correlation.csv")
#Extract network containing significant association (i.e edges) with probability>0.9 (i.e. FDR<0.1)
# net<-extract.network(met.edges, cutoff.ggm = 0.9) 


#--------------------------------------------------------------------------------------------------------
## plot for plasma

library(dplyr)
library(tidyr)
library(ggplot2)



diameter.1 <- 8
nodes.1.P <- read.csv("") %>%  
  filter(Part != "Inorganic") %>% 
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, 5 * pi / 2 - (2 * pi / nrow(.)), length.out = nrow(.))) %>%
  mutate(x = (diameter.1/2) * cos(angle), 
         y = (diameter.1/2) * sin(angle))



diameter.2 <- 2  
nodes.2.P <- read.csv("") %>%
  filter(Part == "Inorganic") %>%
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, 5 * pi / 2 - (2 * pi / nrow(.)), length.out = nrow(.))) %>%
  mutate(x = (diameter.2/2) * cos(angle) + 7, 
         y = (diameter.2/2) * sin(angle) + 2.25)


data.node.P.plot <- rbind(nodes.1.P, nodes.2.P)[-1] 


edges.P <- read.csv("", row.names = 1) 

edges.P <- edges.P %>%
  mutate(From.x = data.node.P.plot$x[match(edges.P$source, data.node.P.plot$label)],
         From.y = data.node.P.plot$y[match(edges.P$source, data.node.P.plot$label)],
         To.x = data.node.P.plot$x[match(edges.P$target, data.node.P.plot$label)],
         To.y = data.node.P.plot$y[match(edges.P$target, data.node.P.plot$label)]) 


#connected_nodes <- unique(c(edges.P$source, edges.P$target))


#data.node.P.plot.connected <- data.node.P.plot[data.node.P.plot$label %in% connected_nodes, ]


ggplot(data = data.node.P.plot) +
  geom_segment(data = edges.P, 
               aes(xend= To.x, yend = To.y, x = From.x, y = From.y, color = factor(correlation), linewidth =factor(correlation)),
               show.legend = T)+
  scale_linewidth_discrete(range = c(0.07, 0.25))+
  geom_point(aes(x = x, y = y, fill = category),
             show.legend = T, size=4, shape = 21, stroke = 0)+ 
  coord_fixed()+
  scale_color_manual(values = c("Negative" = "#427AB2", "Positive" = "#DBDB8D")) +
  scale_fill_manual(values = c("exogenous compound" = "#b5b5d8",
                               "endogenous compound" = "#f6f2c5", 
                               "NEMs" = "#ffcca6", 
                               "EMs" = "#dbcfe4"))+  
  #geom_text(data = data.frame(part = c("Organic", "Inorganic"),
                             # x = c(0,7),
                            #  y = c(0,2.25)),
           # aes(label = part, x = x, y = y), size = 3)+ 
  ## geom_text(data = data.frame(label = c(paste("Nodes:", nrow(data.node.P.plot)),
     ## paste("Edges:", nrow(edges.P)) 
    ##x = 4,
   ## y = c(-3,-3.5,-4)),
    ##aes(x = x, y = y, label = label),hjust = 0, size = 4.1)+
  theme_void() 


ggsave("LEGEND.PNG", width = 5.5, height = 3)

ggsave("Plasma_PLOT2.pdf", width = 5.5, height = 3)

ggsave("Plasma_PLOT2.png", width = 5.5, height = 3)





