if (!require(optparse)) {
  install.packages("optparse")
  library(optparse)
}
if (!require(optparse)) {
  install.packages("tidyverse")
  library(optparse)
}
if (!require(optparse)) {
  install.packages("ggplot2")
  library(optparse)
}
library(ggplot2)
library(tidyverse)
library(optparse)
label_size=10
option_list <- list(
  make_option("--CPMinput", type = "character", default = "/scratch/users/yliu5/HELA_eletter/zenodo/ribo_file/", 
              help = "all analysis folder", metavar = "character"),
  make_option("--GSMinput", type = "character", default = "/scratch/users/yliu5/HELA_eletter/zenodo/HELA_list.csv", 
              help = "input experiment files in csv", metavar = "character"),
  make_option("--outdir", type = "character", default = "/scratch/users/yliu5/HELA_eletter/zenodo/processed/", 
              help = "output directory", metavar = "character")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
CMP_file = "ribo_hela_cpm.csv"
Hela_list= "HELA_list.csv"

Hela_cpm <- read.csv(file.path(args$CPMinput, CMP_file),row.names=1)
infor <- read.csv(file.path(args$GSMinput, Hela_list),row.names=1)[,c("experiment_alias","study_alias")]
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
output_file_path = file.path(args$outdir, "ribo_HeLa_Spearman.csv")
write.csv(hela_table,output_file_path)
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
output_plot_path = file.path(args$outdir, "dedup_ribo_HeLa_spearman_cor.pdf")
ggsave(output_plot_path, width = 2.5, height = 3, units = "in")



