setwd("/Users/yliu5/Library/CloudStorage/Box-Box/Canlab/HELA_science/Hela_fin_data/processed/")
Hela_cpm <- read.csv("ribo_only_cpm_dummy_70.csv",row.names=1)
infor <- read.csv("../HELA_list.csv",row.names=1)[,c("experiment_alias","study_alias")]
out <- data.frame(NULL)
for (i in 1:ncol(Hela_cpm)){
  for (j in i:ncol(Hela_cpm)){
    #print(j)
    m <- cor(Hela_cpm[,i], Hela_cpm[,j],method="spearman")
    out[i,j] <-m
  }}
row.names(out)=colnames(Hela_cpm)
colnames(out)=colnames(Hela_cpm)

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
Hela_cor_tmp <- flattenCorrMatrix(out)
row_tmp <- merge(Hela_cor_tmp,infor,by.x="row",by.y="experiment_alias")
colnames(row_tmp) <- c("row","column","cor","row_id")
col_tmp <- merge(row_tmp,infor,by.x="column",by.y="experiment_alias")
colnames(col_tmp) <- c("row","column","cor","row_id","col_id")
hela_table <- col_tmp %>% mutate(status=ifelse(row_id==col_id,"same_study","dif_study"))
write.csv(hela_table,"hela_table_ribo_only_Spearman.csv")
anno <- hela_table %>% group_by(status) %>% summarise(median_value=median(cor))
ggplot(hela_table,aes(x=status,y=cor,fill=status))+geom_boxplot(alpha=0.8,width=0.3,
                                                                outlier.size = 0.5,outlier.colour = "grey") +
  geom_text(data = anno, aes(y = median_value, label = round(median_value, 3)), 
            vjust = -0.5, color = "black") + 
  theme_bw() +scale_fill_manual(values=c("#984ea3","#4daf4a"))+
  scale_fill_manual(values=c("#984ea3","#4daf4a"))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  theme(axis.text=element_text(size=label_size, colour = "black"),
        axis.title = element_text(size=label_size, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size))+labs(y= "Spearman correlation",
                                                          x="")
ggsave("dedup_ribo_HeLa_spearman_cor.pdf", width = 2.5, height = 3, units = "in")



