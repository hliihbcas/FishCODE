library(ggplot2)
library(ggtree)
args <- commandArgs(trailingOnly = TRUE)
tree = read.tree(args[1])
data=fortify(tree)
tregraph=ggtree(tree, layout="rectangular", size=0.8, col="deepskyblue3") +
geom_tiplab(size=2, color="purple4") + #显示物种信息，并设置颜色大小
#geom_text(aes(subset=!isTip,label=label),size=2, color="deepskyblue3")+#显示所有节点名称
geom_tippoint(size=1, color="deepskyblue3") + #显示物种标识，并设置颜色大小
geom_nodelab(aes(subset=!isTip,label=label),hjust=-0.3, size=2, color="deepskyblue4")+
#geom_text2(aes(subset=!isTip, label=support), hjust=1.3, size=2, color="deepskyblue4") +#显示节点支持率，并设置其位置、大小以及颜色
geom_nodepoint(color="orange", alpha=1/4, size=4) + #显示节点标识及其颜色大小，alpha值为透明度
theme_tree2()+#显示坐标轴（绝对遗传距离）
xlim(NA, max(data$x)*1.2) #调节x轴范围，使得物种信息不超出边界
file_save = paste0(args[1],'.png')
ggsave(file_save,width = 12, height = 8,plot=tregraph)
