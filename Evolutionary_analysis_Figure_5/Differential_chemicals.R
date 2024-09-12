## fold change and FDR

library(tidyverse)
library(rstatix)
rm(list = ls())

data <- read.csv("")
data <- data[-1,]



group <- read.csv("")
group<- group[,-3]
group2<-read.csv("")
group2<- group2[,-3]

## Convert wide data to long format
data1 <- data%>%
  pivot_longer(-X,names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') 

data_group2 <- data%>%
  pivot_longer(-X,names_to = 'Sample', values_to = 'value') %>% 
  left_join(group2, by = 'Sample') 
## FC
data1$value<-as.numeric(data1$value)
dat_FC1 <- data1 %>% group_by(X,Group) %>% 
  summarise(mean=mean(value,na.rm=T),.groups = "keep") %>% 
  pivot_wider(names_from = Group, values_from = mean) %>% 
  summarise("Foldchange_premalignant_IE"=IE/EPL)

data_group2$value<-as.numeric(data_group2$value)
dat_FC2 <- data_group2 %>% group_by(X,Group) %>% 
  summarise(mean=mean(value,na.rm=T),.groups = "keep") %>% 
  pivot_wider(names_from = Group, values_from = mean) %>% 
  summarise("Foldchange_preinvasive_IE"=IE/Pre_invasive)

## ----------------------------------------------------------------------
# 1. EPL VS IE

data[c(-1)]<-lapply(data[c(-1)], as.numeric)
  
## Test for normal distribution

  rownames(data) <- NULL
  dat3 <- data %>% column_to_rownames(var = "X") %>% t() %>% as.data.frame()
  
Vars <- colnames(dat3)
dat_normal_condition <- dat3 %>% 
  shapiro_test(Vars) %>%
  mutate(Normal_Distribution = case_when(p > 0.05 ~ 'YES', TRUE ~ 'NO')) %>% 
  select(X=variable,Normal_Distribution) %>% merge(data,by="X")
head(dat_normal_condition)


dat_wtest1 <- dat_normal_condition %>% 
  filter(Normal_Distribution=="YES") %>% #normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #t_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="EPL",group2=="IE") %>%# A vs B----------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  


dat_wttest2 <- dat_normal_condition %>% 
  filter(Normal_Distribution=="NO") %>% #non-normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #wilcox_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="EPL",group2=="IE") %>%# A vs B------------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  



dat_result <- dat_wtest1 %>% full_join(dat_wttest2) %>% 
  merge(dat_FC1,by="X") %>%
  merge(data,by="X") %>% 
  write_excel_csv(file = "IE_EPLresult.csv")

  
## -------------------------------------------------------------------------
 
# 2.AE+EEVS.IE

dat<-read.csv("", row.names = 1,check.names = F)
df999=reshape2::melt(dat, 
                       id.vars=1, 
                       measure.vars=2:316) 


dat_wtest4 <- df999 %>% group_by(variable) %>% wilcox_test(value ~ label, paired =FALSE ) %>% adjust_pvalue(method = "BH") %>% 
           add_significance("p.adj") %>% filter(group1=="IE",group2=="Pre_invasive") %>% 
          select(-(group1:group2)) %>% 
          select(c('variable','p','p.adj','p.adj.signif')) %>%
          rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif) 

dat_wtest4 <- dat_wtest4 %>%
  mutate(ID = row_number())

data <-data %>%
  mutate(ID = row_number())

dat_result2 <-dat_wtest4 %>%  merge(data,by="ID")
dat_result2 <- dat_result2 [,-2]%>%
  merge(dat_FC2,by="X")%>% 
  write_excel_csv(file = "IE_Pre_invasive_result.csv")

## 3. EE VS IE

data <- read.csv("")
data <- data[-1,]
group <- read.csv("")
group<- group[,-3]
data1 <- data%>%
  pivot_longer(-X,names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') 


group <- read.csv("")
group<- group[,-3]
data1$value<-as.numeric(data1$value)
dat_FC3 <- data1 %>% group_by(X,Group) %>% 
  summarise(mean=mean(value,na.rm=T),.groups = "keep") %>% 
  pivot_wider(names_from = Group, values_from = mean) %>% 
  summarise("Foldchange_EE_IE"=IE/EE)

data[c(-1)]<-lapply(data[c(-1)], as.numeric)

## normal distribution? 

rownames(data) <- NULL
dat4 <- data %>% column_to_rownames(var = "X") %>% t() %>% as.data.frame()

Vars <- colnames(dat4)
dat_normal_condition4 <- dat4 %>% 
  shapiro_test(Vars) %>%
  mutate(Normal_Distribution = case_when(p > 0.05 ~ 'YES', TRUE ~ 'NO')) %>% 
  select(X=variable,Normal_Distribution) %>% merge(data,by="X")
head(dat_normal_condition4)


dat_wtest_EE1 <- dat_normal_condition4 %>% 
  filter(Normal_Distribution=="YES") %>% #normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #t_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="EE",group2=="IE") %>%# A vs B----------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  


dat_wtest_EE2 <- dat_normal_condition4 %>% 
  filter(Normal_Distribution=="NO") %>% #non-normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #wilcox_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="EE",group2=="IE") %>%# A vs B------------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  



dat_result <- dat_wtest_EE1 %>% full_join(dat_wtest_EE2) %>% 
  merge(dat_FC3,by="X") %>%
  merge(data,by="X") %>% 
  write_excel_csv(file = "IE_EEresult.csv")


## 4. health .VS IE

data <- read.csv("")
data <- data[-1,]
group <- read.csv("")
group<- group[,-3]
data1 <- data%>%
  pivot_longer(-X,names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') 



data1$value<-as.numeric(data1$value)
dat_FC4 <- data1 %>% group_by(X,Group) %>% 
  summarise(mean=mean(value,na.rm=T),.groups = "keep") %>% 
  pivot_wider(names_from = Group, values_from = mean) %>% 
  summarise("Foldchange_Hlthy_IE"=IE/Hlthy)

data[c(-1)]<-lapply(data[c(-1)], as.numeric)


rownames(data) <- NULL
dat5 <- data %>% column_to_rownames(var = "X") %>% t() %>% as.data.frame()

Vars <- colnames(dat5)
dat_normal_condition5 <- dat5 %>% 
  shapiro_test(Vars) %>%
  mutate(Normal_Distribution = case_when(p > 0.05 ~ 'YES', TRUE ~ 'NO')) %>% 
  select(X=variable,Normal_Distribution) %>% merge(data,by="X")
head(dat_normal_condition5)


dat_wtest_healthy1 <- dat_normal_condition5 %>% 
  filter(Normal_Distribution=="YES") %>% #normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #t_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="Hlthy",group2=="IE") %>%# A vs B----------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  


dat_wtest_healthy2 <- dat_normal_condition5 %>% 
  filter(Normal_Distribution=="NO") %>% #non-normal distribution
  select(-Normal_Distribution) %>% 
  pivot_longer(-X, names_to = 'Sample', values_to = 'value') %>% 
  left_join(group, by = 'Sample') %>% 
  group_by(X) %>%
  wilcox_test(value ~ Group) %>% #wilcox_test
  adjust_pvalue(method = 'fdr') %>% 
  filter(group1=="Hlthy",group2=="IE") %>%# A vs B------------------------
select(-(group1:group2)) %>% 
  select(c('X','p','p.adj','p.adj.signif')) %>%
  rename(`P_value` = p, FDR = p.adj,signif=p.adj.signif)  



dat_result_healthy <- dat_wtest_healthy1 %>% full_join(dat_wtest_healthy2) %>% 
  merge(dat_FC4,by="X") %>%
  merge(data,by="X") %>% 
  write_excel_csv(file = "IE_healthyresult.csv")
##-----------------------------------------------------------------------
#VOLCANO PLOT

library(ggplot2)
library(ggrepel)

# 1.IE VS epl
data <- read.csv("")



data$Sig <- as.factor(ifelse(data$FDR < 0.05 & data$Foldchange_premalignant_IE > 1, "Up",
                             ifelse(data$FDR < 0.05 & data$Foldchange_premalignant_IE < 1, "Down", "Non")))


write.csv(data,file="")

p2= ggplot(data = data,
           aes(x = Foldchange_premalignant_IE,y = -log10(FDR),colour = Sig,shape= group, fill = Sig))+
  scale_color_manual(values =c("#80B1D3","#bfc0c1","#e88182"))+
  scale_shape_manual(values = c(15, 17, 19,18))+
  geom_point(alpha = 0.8, size = 2)+ 
  scale_size_continuous(range = c(1, 3))+ #指定散点大小渐变
  geom_hline(yintercept =  0 , linetype = "dashed",color= "gray")+
  geom_vline(xintercept = 1, linetype = "dashed",color= "gray")+ 
  theme_bw()+ theme(panel.grid = element_blank())+ #改变主题
  geom_text_repel(aes(label = ""),size = 2.5, max.overlaps = 1000, key_glyph = draw_key_point)+
  labs(x = "Foldchange_premalignant_IE",y = "-log10(FDR)")+
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),# 坐标轴标签和标题
        legend.text = element_text(size = 11),legend.title = element_text(size = 13))+ # 图例标签和标题
  scale_x_continuous(breaks = c(-1,  1, 4), limits = c(-1.7, 4.8))
p2

##-----------------------------------------------------------------------
# 2. IE VS EE+AE

data2 <- read.csv("")


data2$Sig <- as.factor(ifelse(data2$FDR < 0.05 & data2$Foldchange_preinvasive_IE > 1, "Up",
                              ifelse(data2$FDR < 0.05 & data2$Foldchange_preinvasive_IE < 1, "Down", "Non")))
write.csv(data2,file="")

p3= ggplot(data = data2,
           aes(x = Foldchange_preinvasive_IE,y = -log10(FDR),colour = Sig,fill = Sig,shape= group))+
  scale_color_manual(values =c("#80B1D3","#bfc0c1","#e88182"))+
  scale_shape_manual(values = c(15, 17, 19,18))+
  geom_point(alpha = 0.8, size = 2)+ 
  scale_size_continuous(range = c(1, 3))+ 
  geom_hline(yintercept =  0 , linetype = "dashed",color= "gray")+
  geom_vline(xintercept = 1, linetype = "dashed",color= "gray")+ 
  theme_bw()+ theme(panel.grid = element_blank())+ 
  geom_text_repel(aes(label = ""),size = 2.5, max.overlaps = 1000, key_glyph = draw_key_point)+
  labs(x = "Foldchange_preinvasive_IE",y = "-log10(FDR)")+
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),# 坐标轴标签和标题
        legend.text = element_text(size = 11),legend.title = element_text(size = 13))+ 
  xlim(c(-0.5, 3))

p3

##-----------------------------------------------------------------------------------

# 3.IE VS EE

data3 <- read.csv("")



data3$Sig <- as.factor(ifelse(data3$FDR < 0.05 & data3$Foldchange_EE_IE > 1, "Up",
                              ifelse(data3$FDR < 0.05 & data3$Foldchange_EE_IE < 1, "Down", "Non")))

ggplot(data = data3,
       aes(x = Foldchange_EE_IE,y = -log10(FDR),colour = Sig,fill = Sig,shape= group))+
  scale_color_manual(values =c("#80B1D3","#bfc0c1","#e88182"))+
  scale_shape_manual(values = c(15, 17, 19,18))+
  geom_point(alpha = 0.8, size = 2)+
  scale_size_continuous(range = c(1, 3))+ 
  geom_hline(yintercept =  0 , linetype = "dashed",color= "gray")+
  geom_vline(xintercept = 1, linetype = "dashed",color= "gray")+ 
  theme_bw()+ theme(panel.grid = element_blank())+ 
  geom_text_repel(aes(label = ""),size = 2.5, max.overlaps = 1000, key_glyph = draw_key_point)+
  labs(x = "Foldchange_EE_IE",y = "-log10(FDR)")+
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13))+ 
  xlim(c(-0.5, 3))


# 4.IE VS HC

data4 <- read.csv("")



data4$Sig <- as.factor(ifelse(data4$FDR < 0.05 & data4$Foldchange_Hlthy_IE > 1, "Up",
                              ifelse(data4$FDR < 0.05 & data4$Foldchange_Hlthy_IE < 1, "Down", "Non")))
write.csv(data4,file="IE_healthresult_final.csv")

p1= ggplot(data = data4,
           aes(x = log10(Foldchange_Hlthy_IE),y = -log10(FDR),colour = Sig,fill = Sig,shape= group))+
  scale_color_manual(values =c("#80B1D3","#bfc0c1","#e88182"))+
  scale_shape_manual(values = c(15, 17, 19,18))+
  geom_point(alpha = 0.8, size = 2)+ 
  scale_size_continuous(range = c(1, 3))+ 
  geom_hline(yintercept =  0 , linetype = "dashed",color= "gray")+
  geom_vline(xintercept = 0, linetype = "dashed",color= "gray")+ 
  theme_bw()+ theme(panel.grid = element_blank())+ 
  geom_text_repel(aes(label = ""),size = 2.5, max.overlaps = 1000, key_glyph = draw_key_point)+
  labs(x = "log10_Foldchange_Hlthy_IE",y = "-log10(FDR)")+
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),# 坐标轴标签和标题
        legend.text = element_text(size = 11),legend.title = element_text(size = 13))
p1


library(patchwork)
p1_no_legend <- p1 + guides(color = "none", shape = "none", fill = "none")
p2_no_legend <- p2 + guides(color = "none", shape = "none", fill = "none")


combined_plot <- p1_no_legend + p2_no_legend +  p3


print(combined_plot)

ggsave(filename = "volcano.pdf", plot = combined_plot, width = 11, height = 3.5)



##-------------------------------------------------------------------------------------------
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(Seurat)
library(ggprism)
library(patchwork)

rm(list= ls() ) 
## -----------------------------------------------------------------------------------------------------------------
## Kw test with BH correction

library(stats)

T_Chemicals <- read.csv("", check.names = FALSE)

# Create an empty list for storing Kruskal-Wallis and BH correction results
results_list <- list()

# Create a vector to store variable names
variable_names <- colnames(T_Chemicals)[3:317]

for (i in 3:317) {
  # Execute Kruskal-Wallis test
  kw_result <- kruskal.test(T_Chemicals[, i], as.factor(T_Chemicals[, 2]))
  
  # Store Kruskal-Wallis test result in the list along with variable name
  result_entry <- list(Variable = variable_names[i - 2], Kruskal_Test = kw_result)
  results_list[[i - 2]] <- result_entry
}

# Extract p-values from the Kruskal-Wallis test results
p_values <- sapply(results_list, function(result) result$Kruskal_Test$p.value)

# Apply BH correction to the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Add the adjusted p-values to the results list
for (i in 1:length(results_list)) {
  results_list[[i]]$Kruskal_Test$adjusted_p_value <- adjusted_p_values[i]
}

# Write the results to a text file
capture.output(results_list, file = "kw_with_BH.list.txt", append = FALSE)

##--------------------------------------------------------------------------------------------------
## plot all significant altered elements
df_chemicals = reshape2::melt(T_Chemicals, Id.vars=2, measure.vars=3:317)



df_chemicals = read.csv("")



my_theme <- theme_prism() +theme(axis.line = element_line(size = 1),
                                 panel.border = element_rect(colour = "black", fill = NA, size = 1),
                                 axis.text.x = element_text(size = 8, angle = 40))



ggplot(data = df_chemicals, aes(x = factor(Stage,levels=c("hlthy", "EPL","Tis/T1", "T2", "T3", "T4")), y = value, fill = Stage))+ ##用stage列设置散点图和箱线图的填充颜色
  geom_violin(position = position_dodge(0.8),trim = T)+
  stat_summary(fun = mean, geom = "point", shape = 17, size = 0.8, color = "black")+
  facet_wrap( ~ variable,nrow=3, scales = "free") +
  my_theme+ scale_fill_manual(values = c("#1c9e77","#D26900","#FFFACD","#467500","#fdae6b", "#FFDAB9"))

ggsave("", width = 10, height = 10)



ggplot(data = df_chemicals, aes(x = factor(Stage,levels=c("hlthy", "EPL","Tis/T1", "T2", "T3", "T4")), y = value, fill = Stage))+ ##用stage列设置散点图和箱线图的填充颜色
  geom_violin(position = position_dodge(0.8),trim = T)+
  stat_summary(fun = mean, geom = "point", shape = 17, size = 0.8, color = "black")+
  scale_y_log10() +
  facet_wrap( ~ variable,nrow=3, scales = "free") +
  my_theme+ scale_fill_manual(values = c("#1c9e77","#D26900","#FFFACD","#467500","#fdae6b", "#FFDAB9"))

ggsave("", width = 10, height = 10)

##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

## RIDGE PLOT

library(ggridges)
library(ggplot2)
library(cols4all)

ridge_meta = read.csv("",check.names = F)

df_ridge_meta = reshape2::melt(ridge_meta, Id.vars=1, measure.vars=2:300)
y_axis_order <- c("hlthy", "EPL","Tis/T1", "T2", "T3", "T4")
mycol <- c4a('superfishel_stone',6)

p <- ggplot(data = df_ridge_meta,
            aes(x = value, y = stage, fill = stage)) +
  geom_density_ridges(alpha = 0.8,color = 'white') +
  theme_classic() +
  scale_y_discrete(limits = y_axis_order)+
  scale_fill_manual(values = mycol)

p

##--------------------------------------------

ridge_meta2 = read.csv("",check.names = F)

df_ridge_meta2 = reshape2::melt(ridge_meta2, Id.vars=1, measure.vars=2:228)
y_axis_order <- c("hlthy", "EPL","Tis/T1", "T2", "T3", "T4")
mycol <- c4a('superfishel_stone',6)

p2 <- ggplot(data = df_ridge_meta2,
             aes(x = value, y = stage, fill = stage)) +
  geom_density_ridges(alpha = 0.8,color = 'white') +
  theme_classic() +
  scale_y_discrete(limits = y_axis_order)+
  scale_fill_manual(values = mycol)

p2
##--------------------------------------------

ridge_meta3 = read.csv("",check.names = F)

df_ridge_meta3 = reshape2::melt(ridge_meta3, Id.vars=1, measure.vars=2:73)
y_axis_order <- c("hlthy", "EPL","Tis/T1", "T2", "T3", "T4")
mycol <- c4a('superfishel_stone',6)

p3 <- ggplot(data = df_ridge_meta3,
             aes(x = value, y = stage, fill = stage)) +
  geom_density_ridges(alpha = 0.8,color = 'white') +
  theme_classic() +
  scale_y_discrete(limits = y_axis_order)+
  scale_fill_manual(values = mycol)

p3


