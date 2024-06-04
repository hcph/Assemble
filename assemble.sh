##确定contig的顺序，将contig利用minimap2比对到参考基因组上，每个材料比对一遍
bsub  -J asm -n 5 -o asm.out -e asm.err -q smp -R "rusage[mem=300GB]" "minimap2 -t 5 -x asm5 2DL.fa YZ1_2dl.p_ctg.fa > YZ1_2dl_totalsample.paf"
#YZ1
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+ArinaLrFor.fa  YZ1_2dl.p_ctg.fa > YZ1_ArinaLrFor.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+CDC_Landmark.fa YZ1_2dl.p_ctg.fa > YZ1_CDC_Landmark.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+CDC_Stanley.fa YZ1_2dl.p_ctg.fa > YZ1_CDC_Stanley.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+IWGSC.fa YZ1_2dl.p_ctg.fa > YZ1_IWGSC.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Jagger.fa YZ1_2dl.p_ctg.fa > YZ1_Jagger.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Julius.fa YZ1_2dl.p_ctg.fa > YZ1_Julius.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Lancer.fa YZ1_2dl.p_ctg.fa > YZ1_Lancer.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Mace.fa YZ1_2dl.p_ctg.fa > YZ1_Mace.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Norin.fa YZ1_2dl.p_ctg.fa > YZ1_Norin.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+PL190962.fa YZ1_2dl.p_ctg.fa > YZ1_PL190962.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_AK58.fa YZ1_2dl.p_ctg.fa > YZ1_AK58.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_KN9204.fa YZ1_2dl.p_ctg.fa > YZ1_KN9204.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_fielder.fa YZ1_2dl.p_ctg.fa > YZ1_fielder.paf
#YM33
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+ArinaLrFor.fa  YM33_2dl.p_ctg.fa > YM33_ArinaLrFor.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+CDC_Landmark.fa YM33_2dl.p_ctg.fa > YM33_CDC_Landmark.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+CDC_Stanley.fa YM33_2dl.p_ctg.fa > YM33_CDC_Stanley.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+IWGSC.fa YM33_2dl.p_ctg.fa > YM33_IWGSC.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Jagger.fa YM33_2dl.p_ctg.fa > YM33_Jagger.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Julius.fa YM33_2dl.p_ctg.fa > YM33_Julius.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Lancer.fa YM33_2dl.p_ctg.fa > YM33_Lancer.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Mace.fa YM33_2dl.p_ctg.fa > YM33_Mace.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+Norin.fa YM33_2dl.p_ctg.fa > YM33_Norin.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_10+PL190962.fa YM33_2dl.p_ctg.fa > YM33_PL190962.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_AK58.fa YM33_2dl.p_ctg.fa > YM33_AK58.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_KN9204.fa YM33_2dl.p_ctg.fa > YM33_KN9204.paf
minimap2 -x asm5  /public/home/ynshen/pangenome/modified_chr2D_fielder.fa YM33_2dl.p_ctg.fa > YM33_fielder.paf
##提取mapq = 60的contig比对上的位置,转化为bed文件
cut -f 1 /public/home/tllu/2DL_HiFi/pan_blast/2DL/sample.txt | while read i;
do
awk '{if($12 >= 60) print "chr2D\t"$8"\t"$9"\t"$1}'  YM33_"$i".paf | sort -k1,1 -k2n,2 > visual/YM33_MAPQ60_"$i".bed
#合并bed文件
cut -f 1 /public/home/tllu/2DL_HiFi/pan_blast/2DL/sample.txt | while read i;
do
bedtools merge -i YM33_MAPQ60_"$i".bed -c 4 -o collapse > YM33_MAPQ60_"$i".merged.bed
done
#统计每个材料的在500000000-600000000的覆盖度
cut -f 1 /public/home/tllu/2DL_HiFi/pan_blast/2DL/sample.txt | while read i;
do
bedtools closest -D ref -t all -mdb all -a YM33_MAPQ"$j"_"$i".merged.bed -b target_region.bed \
|  awk  '{if($NF==0)print $1"\t"$2"\t"$3"\t"$4}'  > YM33_MAPQ"$j"_"$i".merged_2DL.bed
done
#统计scaffold长度
wc -l YM33_MAPQ60_*.merged_2DL.bed >YM33_MAPQ60_scaffold_length.txt
####可视化覆盖度
###R
##绘制不同材料的覆盖度
library(pafr)
library(ggplot2)
library(tidyverse)
library(patchwork)
#导入数据
setwd("D:/2DL三代测序/2DL/YZ1")
name <- c("AK58","ArinaLrFor", "CDC_Landmark", "CDC_Stanley", "IWGSC", "Jagger", "Julius","KN9204","Lancer", "Mace", "Norin", "PL190962","fielder")
## 查看结果
for (j in c(0,10,20,30,40,50,60)) {
  plot_list <- list()
  for (i in name) {
    df <- read_paf(paste0("YZ1_",i,"_MAPQ",j,".paf"))
    p1 <- plot_coverage(df,fill='orange') + 
      theme(legend.position = "none") + 
      xlim(460000000,630000000) + 
      ylab(paste0(">=","MAPQ",j)) +
      xlab(NULL) + 
      ggtitle(i)
    plot_list[[i]] <- p1
  }
  #组图
  grid.arrange(plot_list[[2]],plot_list[[6]], 
               plot_list[[7]],plot_list[[13]], 
               plot_list[[1]],plot_list[[3]],
               plot_list[[4]],plot_list[[5]],
               plot_list[[8]],plot_list[[9]],
               plot_list[[10]],plot_list[[11]],
               plot_list[[12]],
               nrow=13,ncol=1)     %>%  ggsave(paste0("YZ1_MAPQ",j,"_coverage.pdf"),.,width=100,height=300, units="mm")
}
setwd("D:/2DL三代测序/2DL/YM33")
for (j in c(0,10,20,30,40,50,60)) {
  plot_list <- list()
  for (i in name) {
    df <- read_paf(paste0("YM33_",i,"_MAPQ",j,".paf"))
    p1 <- plot_coverage(df,fill='orange') + 
      theme(legend.position = "none") + 
      xlim(460000000,630000000) + 
      ylab(paste0(">=","MAPQ",j)) +
      xlab(NULL) + 
      ggtitle(i)
    plot_list[[i]] <- p1
  }
  #组图
  grid.arrange(plot_list[[2]],plot_list[[6]], 
               plot_list[[7]],plot_list[[13]], 
               plot_list[[1]],plot_list[[3]],
               plot_list[[4]],plot_list[[5]],
               plot_list[[8]],plot_list[[9]],
               plot_list[[10]],plot_list[[11]],
               plot_list[[12]],
               nrow=13,ncol=1)     %>%  ggsave(paste0("YZ1_MAPQ",j,"_coverage.pdf"),.,width=100,height=300, units="mm")
}
