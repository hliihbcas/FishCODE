#!/home/hli/miniconda3/envs/visualomics_02/bin/Rscript
library(ggplot2)
library(rjson)
#./Rversion1.R ordination 11 ./100100100/*.data 0 ./100100100/feature 2> ./100100100/err_11.txt 
#./Rversion1.R 0 01 0 ./100100100/user.conf 0 2> ./100100100/err_01.txt 

main <- function(){
    argstr<- commandArgs(T)
    type <- argstr[1]#top type
	subtype <- argstr[2]#sub type
    datafile <- argstr[3]#不一定能成
    plot.conf <- argstr[4]
	feature.conf <- argstr[5]
	#feature <- fromJSON(file=feature.conf)
    if (subtype == "11"){
		if (type == "de_analyse_batch_time"){
			library(ggraph)
			library(tidyverse)
			library(netET)
			library(RColorBrewer)
			feature <- fromJSON(file=feature.conf)
			#11 重读文件
					
			#write.table(df,scole_file,header=T)
			lev <- read.table(paste0(feature$out,"/Node_level.tsv"), header=T,com='', quote='',check.names=F, sep="\t")
			names(lev)=c('TF.gene.ID','level.in.GCN')
			print(str(lev))
			

			lev2 <- lev %>% mutate(level.in.GCN = factor(paste0("lev",level.in.GCN))) 
			lev2$level.in.GCN <- fct_relevel(lev2$level.in.GCN,paste0("lev",1:30))
			#lev2$assigned_level <- fct_relevel(lev2$assigned_level,paste0("lev",1:16))
			head(lev2)

			TF_gene <- read.table(paste0(feature$out,"/TF_gene.txt"), header=F,com='', quote='',check.names=F, sep="\t",row.names=1)
			
			TF_gene = as.data.frame(t(TF_gene))
			cor_matrix <- cor(TF_gene)
			threshold <- as.numeric(feature$cut_level)
			significant_pairs <- which(cor_matrix >= threshold & upper.tri(cor_matrix, diag = FALSE), arr.ind = TRUE)
			tt = apply(significant_pairs,1,TF_gene_a,cor_matrix)
			pcc = as.data.frame(t(tt))
			colnames(pcc) = c('TF.gene.ID','gene.ID','p_value')

			pcc %>%
    			filter(TF.gene.ID %in% lev$TF.gene.ID, gene.ID %in% lev$TF.gene.ID) %>%
    			slice_head(n=20000) -> pcc2
			
			level_counts <- lev2 %>% 
    			filter(TF.gene.ID %in% unique(c(pcc2$gene.ID ,pcc2$TF.gene.ID))) %>%  
    			count(level.in.GCN) %>% 
    			as.data.frame()

			lev3 <- lev2 %>%
    			left_join(level_counts, by = "level.in.GCN") %>%
    			mutate(label.in.GCN = paste0(level.in.GCN, " (n=", n, ")"))

			#pcc <- read.csv(paste0(feature$out,"/Node_level.tsvNode_relation.csv"), header=T,com='', quote='',check.names=F)
			as_tbl_graph(pcc2[,c("TF.gene.ID","gene.ID")],simplify = "all") %>% left_join(lev3, by = c("name" = "TF.gene.ID"))  ->g

			if (feature$layout=="Circle"){
				coord <- layout_on_circle(g, group = level.in.GCN)
			}
			if (feature$layout=="Horizontal"){
				n_levels = length(unique(lev$level.in.GCN))
				x_range <- seq(from = -1, to = 1, length.out = n_levels)
				y_range <- seq(from = 0, to = 0, length.out = n_levels)

				labels <- paste0("lev", 1:n_levels)

				center <- list(x = setNames(x_range, labels),y = setNames(y_range, labels))
				print(center)
				coord <- layout_on_circle(g, group = level.in.GCN, center = center)
			}
			if (feature$layout=="45_degrees"){
				n_levels = length(unique(lev$level.in.GCN))
				x_range <- seq(from = -1, to = 1, length.out = n_levels)
				y_range <- seq(from = -1, to = 1, length.out = n_levels)

				labels <- paste0("lev", 1:n_levels)

				center <- list(x = setNames(x_range, labels),y = setNames(y_range, labels))
				print(center)
				coord <- layout_on_circle(g, group = level.in.GCN, center = center)
			}


			p <- ggraph(g, coord, x = coord[,1],y = coord[,2]) +

    			geom_edge_fan(width = 0.2, alpha = 0.1) +

    			geom_node_point(aes(colour = level.in.GCN),
                    size = 0.5,
                    show.legend = F) +
                     
    			coord_fixed(clip = "off") +
    
    			theme(panel.background = element_blank()) +

    			ggforce::geom_mark_circle(
        			aes(x,  
            			y, 
            			group = level.in.GCN, 
            			label = label.in.GCN),
        			label.fontsize = 10,
        			con.type = "straight",
        			label.fontface = c("bold", "italic"),
        			color = NA,
        			expand = unit(1, "mm"),
        			alpha = 0.25 ) +
        
    				scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(30))

			
	
			feature$mutiformat='valcano'
			plot_continue(p,feature)

			second_step <- read.table(paste0(feature$out,"/TF_gene.txt"), header=F,com='', quote='',check.names=F, sep="\t")
			print(str(lev))
			print(str(second_step))
			second_step_col = colnames(second_step)[2:length(colnames(second_step))]
			second_step_num = second_step[,colnames(second_step)[2:length(colnames(second_step))]]
			second_step_num_t = as.data.frame(t(apply(second_step_num,1,scale)))
			colnames(second_step_num_t) = second_step_col
			second_step = data.frame(V1=second_step$V1,second_step_num_t)
			print(str(second_step))
			new_heatmap = merge(lev,second_step,by.x = "TF.gene.ID", by.y = "V1")
			print('gg')
			print(str(new_heatmap))
			max_len = max(new_heatmap$level.in.GCN)
			df = list()
			for (i in sort(unique(new_heatmap$level.in.GCN))){
				sub_dataframe = subset(new_heatmap,level.in.GCN==i)
				lev_name = colnames(sub_dataframe)[3:ncol(sub_dataframe)]
				pre_sub = subset(sub_dataframe,select=lev_name)
				result_t = sapply(pre_sub,mean)
				labels = paste0("level",i)
				df[[labels]] = result_t
			}
			print(str(df))
			df_new = as.data.frame(t(as.data.frame(df)))


			conditions_f <- read.table(paste0(feature$out,"/mysql_condition.txt"),com='', quote='',check.names=F, sep="\t")
			pre_colname = sort(unique(conditions_f$condition))
			
			colnamess = c()
			for(i in pre_colname){
				temp_num = paste0("time point " ,strsplit(i,'_')[[1]][2])

				colnamess = c(colnamess,temp_num)
			}

			
			colnames(df_new) = colnamess
			print(str(df_new))
			print(head(df_new))
			feature$mutiformat='pca'

			feature$process = '0'
			feature[['r_method']] = 'pheatmap'
			feature$legend = "yes"
			feature$number = "no"
			feature$title = ""
			feature$angle = c(45)
			feature$size = c(10)

			print(feature)
			heatmap_pheatmap(df_new,feature)

			

			
			
		}
		if (type == "ordinations"){
			feature <- fromJSON(file=feature.conf)
			#11 重读文件
			data <- read.delim(datafile,stringsAsFactors=F,blank.lines.skip = TRUE,comment.char = '')			
			#write.table(df,scole_file,header=T)
			if (feature$method1 == "pcoa"){
				#sre.png,biplot.png,res,res的value和score文件以及包含前俩个主成分的文件可以用于其他的绘图，包括group以及添加了一列每行的名字
				#ordination.png 以及对应的ggplot2的保留文件
				pcoa_main(data,feature,subtype)
			}
		}
		if(type == 'barplots'){
			tryCatch({
			feature <- fromJSON(file=feature.conf)
			data <- read.delim(datafile,stringsAsFactors=F)
			data[[feature$group]]=factor(data[[feature$group]],levels=unique(data[[feature$group]]))
			
			if(feature$bar_type == 'mean'){
				#data[[feature$group]]=factor(data[[feature$group]],levels=unique(data[[feature$group]]))
				print(levels(data[[feature$group]]))
				bar_mean(data,feature,subtype)
			}
			else{
data[[feature$subgroup]]=factor(data[[feature$subgroup]],levels=unique(data[[feature$subgroup]]))
				bar_charts(data,feature,subtype)
			}
			},error=function(e){
				print("aaa")
				p <- ggplot(mtcars,aes(gear,carb)) + theme_void() + annotate("text",x=range(mtcars$gear)[1],y=range(mtcars$carb)[1],label="FishCODE detected an operation anomaly. \nPlease contact the author promptly.",size=6)
				print("asdas")
				#feature$sample_name=null
				plot_continue(p,feature)
				
				return("1")
			})

		}
		if(type== 'batch_boxplots'){
			feature <- fromJSON(file=feature.conf)
			#feature$stand = 'yes'
			

			tryCatch({
			feature$stand = 'no'
			feature$group = "group"
			feature$method1 = "pcoa"
			feature$distance = "euclidean"
			feature$typed = "t"
			feature$level = "0.95"
			feature$count=0
			data <- read.table(datafile, header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
			if (file.exists(paste0(feature$out,"/mysql_selected.txt"))){
				dara_02 = read.delim(paste0(feature$out,"/mysql_selected.txt"),stringsAsFactors=F)
			}else{
				#print(paste0("cp ",datafile," ",feature$out,"/mysql_selected.txt"))
				system(paste0("cp ",datafile," ",feature$out,"/mysql_selected.txt"))
				dara_02 = read.table(paste0(feature$out,"/mysql_selected.txt"), header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
				#read.delim(paste0(feature$out,"/mysql_selected.txt"),stringsAsFactors=F)
			}
			
			print(colnames(dara_02));print(colnames(data))
			if(mean(colnames(data) == colnames(dara_02))!=1){
				insert_names = intersect(rownames(data),rownames(dara_02))
				data_expr = cbind(data[insert_names,],dara_02[insert_names,])
			}else{
				data_expr = dara_02
			}
			data_expr = na.omit(data_expr)
			data_expr_save = apply(data_expr,2,round,0)
			save_file_name = paste0(feature$out,"/end_expr.txt")
			
			write.table(data_expr_save,save_file_name,sep='\t',quote=FALSE)


			condition = read.delim(paste0(feature$out,"/mysql_condition.txt"),stringsAsFactors=F)
			#insert_names = intersect(rownames(data),rownames(dara_02))
			#print(insert_names)
			#data_expr = cbind(data[insert_names,],dara_02[insert_names,])
			data_pre = as.data.frame(t(data_expr))
			data_pre$group = condition[['Project.name']]
			print(condition)
			print(str(data_pre));print(rownames(data_pre))
			#write.table(data_mapped[c('STRING_id','SYMBOL')],data_mapped_file,sep='\t',row.names=F,col.names=F,quote=FALSE)

			},error=function(e){
				print("aaa")
				p <- ggplot(mtcars,aes(gear,carb)) + theme_void() + annotate("text",x=range(mtcars$gear)[1],y=range(mtcars$carb)[1],label="Detected abnormal input data and parameter reading,\n please refer to the standard input file for proofreading.\nPlease browse the Err report panel at the bottom of the page, \nwhich contains detailed error information.",size=6)
				feature$sample_name=null
				plot_continue(p,feature)
				
				return("1")
			})
			pcoa_main(data_pre,feature,subtype)

		}
		if(type=='batch_boxplots_batch'){
			

			#tryCatch({
			library("sva")
			feature <- fromJSON(file=feature.conf)
			#feature$stand = 'no'
			feature$group = "group"
			feature$method1 = "pcoa"
			feature$distance = "euclidean"
			feature$typed = "t"
			feature$level = "0.95"
			feature$count=0
			data <- read.table(datafile, header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
			dara_02 = read.delim(paste0(feature$out,"/mysql_selected.txt"),stringsAsFactors=F)
			print(colnames(data))
			print(colnames(dara_02))
			print(mean(colnames(data) != colnames(dara_02)))
			if(mean(colnames(data) != colnames(dara_02))!=0){
				
				insert_names = intersect(rownames(data),rownames(dara_02))
				data_expr = cbind(data[insert_names,],dara_02[insert_names,])
			}else{
				data_expr = dara_02
			}
			condition = read.delim(paste0(feature$out,"/mysql_condition.txt"),stringsAsFactors=F)
			batch  = condition[['Project.name']]
			tissue = condition[['condition']]
			mod = model.matrix(~as.factor(condition),data=condition)
			print(table(batch,tissue))
			
			save_row = rownames(data_expr)
			save_col = colnames(data_expr)
			print(head(data_expr))

			#数据清洗
			data_expr = na.omit(data_expr)
			#data_expr = as.data.frame(lapply(data_expr,round,0))

			#低表达量删选
			if(feature['batch']=="ComBat" || feature['batch']=="DEseq"){
				data_expr = data_expr[rowSums(data_expr) >= length(data_expr),]
			}

			#数据标准化
			data_expr = remove_batch_stand(data_expr,feature)
			
			print(head(data_expr));print(class(data_expr))
			#标准误筛选
			data_expr = data_expr[apply(data_expr,1,var)>=feature['sd'],]
			print(head(data_expr));print(str(data_expr));print(class(data_expr))
			#取整
			data_expr = apply(data_expr,2,round,0)
		
			#save_file_name = paste0(feature$out,"/before_expr.txt")
			
			
			#write.table(data_expr,save_file_name,sep='\t',quote=FALSE)
			
			print(head(data_expr));print(str(data_expr));print(class(data_expr))
			#批次效应删除
			if(feature['batch']=="ComBat"){
				library("BiocParallel")
				register(MulticoreParam(4))
				#data_expr = ComBat_seq(dat=data_expr,batch=batch,mod=mod,par.prior=T, prior.plots=F)
				data_expr = ComBat_seq(data_expr,batch=batch,group=tissue, full_mod = TRUE)
				print(str(data_expr));print(head(data_expr))
			}
			if(feature['batch']=="DEseq"){
				save_file_name = paste0(feature$out,"/end_expr.txt")
				write.table(data_expr,save_file_name,sep='\t',quote=FALSE)
				condition = condition[colnames(data_expr),,drop=F]
				denanlyse_batch(data_expr,condition,feature)
				return('mad')
			}
			
			save_file_name = paste0(feature$out,"/end_expr.txt")
			
			
			write.table(data_expr,save_file_name,sep='\t',quote=FALSE)
			
			data_expr = as.data.frame(t(data_expr))
			data_expr$group = condition[['Project.name']]
			feature$stand = 'no'

			
			#},error=function(e){
			#	p <- ggplot(mtcars,aes(gear,carb)) + theme_void() + annotate("text",x=range(mtcars$gear)[1],y=range(mtcars$carb)[1],label="Detected abnormal input data and parameter reading,\n please refer to the standard input file for proofreading.\nPlease browse the Err report panel at the bottom of the page, \nwhich contains detailed error information.",size=6)
			#	feature$sample_name=null
			#	plot_continue(p,feature)
			#	return("1")
			#})



			pcoa_main(data_expr,feature,subtype)
		}
		if(type == 'boxplots'){
			feature <- fromJSON(file=feature.conf)
			data <- read.delim(datafile,stringsAsFactors=F)
			data[[feature$group]]=factor(data[[feature$group]],levels=unique(data[[feature$group]]))			
			if(feature$bar_type == 'mean'){
				box_mean(data,feature,subtype)
			}
			else{
				data[[feature$subgroup]]=factor(data[[feature$subgroup]],levels=unique(data[[feature$subgroup]]))
				box_charts(data,feature,subtype)
			}
		}

		if(type == 'heatmap'){
			feature <- fromJSON(file=feature.conf)
			data <- read.delim(datafile,stringsAsFactors=F,row.names = 1)
			if(feature[['r_method']]=='ggplot2'){
				heatmap_all(data,feature)
			}
			else{
				heatmap_pheatmap(data,feature)
			}
		}
		if(type == 'de_analyse'){
			feature <- fromJSON(file=feature.conf)
			data <- read.table(datafile, header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
			table_path = paste0(feature$out,'/sample')
			table <- read.table(table_path,header = TRUE,sep = "\t",row.name=1)
			#data <- read.delim(datafile,stringsAsFactors=F,row.names = 1)
			deanalyse(data,table,feature)
		}
		if(type == 'de_analyse_batch'){
			print('aa')
			feature <- fromJSON(file=feature.conf)
			

			tryCatch({
				print('aaa')
				print(paste0(feature$out,"/end_expr.txt"))
				data <- read.table(paste0(feature$out,"/end_expr.txt"), header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
				print('a')
			table_path = paste0(feature$out,'/mysql_condition.txt')
			table <- read.delim(table_path)
			table[,2] = factor(table[,2],levels=unique(table[,2]))
			#table[['Project.name']] = NULL

			data = data[,rownames(table),drop=F]
			print(str(data))
			print(str(table))
			#table_new = table[colnames(data),,drop=F]
			table_new = table
			#data <- read.delim(datafile,stringsAsFactors=F,row.names = 1)
			
			#data = as.data.frame(lapply(data,round,0))
			#data = na.omit(data)
			print("b");print(str(data));print(head(data))
			},error=function(e){
				p <- ggplot(mtcars,aes(gear,carb)) + theme_void() + annotate("text",x=range(mtcars$gear)[1],y=range(mtcars$carb)[1],label="Detected abnormal input data and parameter reading,\n please refer to the standard input file for proofreading.\nPlease check whether you upload the necessary group files\nIf it still does not work properly after checking, \nPlease browse the Err report panel at the bottom of the page, \nwhich contains detailed error information.",size=6)
				feature$sample_name=null
				feature$mutiformat='heat'
				feature$r_method = 'ggplot2'
				plot_continue(p,feature)

				feature$mutiformat='valcano'
				plot_continue(p,feature)

				feature$mutiformat='pca'
				plot_continue(p,feature)

				return("1")
			})



			deanalyse(data,table_new,feature)
		}

		if(type == 'enrich_analyse'){
			feature <- fromJSON(file=feature.conf)
			data <- read.delim(datafile,stringsAsFactors=F)
			enrich_analyse(data,feature)
		}
		if(type == 'ppi'){
			feature <- fromJSON(file=feature.conf)
			data <- read.delim(datafile,stringsAsFactors=F)
			ppi_analyse(data,feature)
		}

		if(type == 'domain'){
			#library(gggenes)
			library(drawProteins)
			feature <- fromJSON(file=feature.conf)
			input_sub_data = read.table(datafile,head=T,stringsAsFactors=F,sep='\t')
			#input_data = subset(input_sub_data,select=c('molecule','gene','start','end','strand','orientation'),drop=F)
			#input_data = input_data[!duplicated(input_data),]
			#print(str(input_data))
			#print(str(input_sub_data))
			#a = apply(input_data,1,gene_plot,out=feature,gene=input_data,sub_gene=input_sub_data)			
			domain_gene_plot(input_sub_data,feature)
		}

	}

	if (subtype == "01"){
		detail_plot(plot.conf)
	}

}


TF_gene_a = function(x,cor_matrix){
	source_name = colnames(cor_matrix)[x[1]]
	target_name = colnames(cor_matrix)[x[2]]
	pcc_value = cor_matrix[x[1],x[2]]
	if(pcc_value)
	
	return(c(source_name,target_name,pcc_value))
}

minmax = function(x){
	minx = min(x)
	maxx = max(x)
	return_x = c()
	for(i in x){	
		return_x=c(return_x,(i-minx)/(maxx-minx))
	}
	return(return_x)

}

meanstand = function(x){
	meanx = mean(x)
	maxx = max(x)
	minx = min(x)
	return_x = c()
	for(i in x){
		return_x = c(return_x,(i-meanx)/(maxx-minx))
	}
	return(return_x)
}

remove_batch_stand <- function(data,feature){
	print(feature['stand'])
	col_name = colnames(data)
	if (feature['stand']=="no"){
		return(data)
	}
	if(feature['stand']=="log"){
		data = as.data.frame(t(apply(data+1,1,log)))
		
		return(data)
	}
	if(feature['stand']=="minmax"){
		data = as.data.frame(t(apply(data,1,minmax)))
		colnames(data) = col_name
		return(data)
	}
	if(feature['stand']=="mean"){
		data = as.data.frame(t(apply(data,1,meanstand)))
		colnames(data) = col_name
		return(data)
	}
	if(feature['stand']=="standard"){
		data = as.data.frame(t(apply(data,1,scale,center=T,scale=T)))
		colnames(data) = col_name
		return(data)
	}


	
}


domain_gene_plot <- function(data,feature){
	p <- draw_canvas(data)
	p<- draw_chains(p,data)
	#p <- draw_domains(p,data)
	draw_type = unique(data$type)
    for(i in draw_type){
		if(i !='CHAIN'){
			p <- draw_domains(p,data,type=i)
		}
	}
	p <- p + theme_bw(base_size = 20) + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +theme(axis.ticks =element_blank(),axis.text.y =element_blank()) + theme(panel.border =element_blank())
	plot_continue(p,feature)
	
}

gene_plot <- function(x,out,gene,sub_gene){

    name = strsplit(x[1],':')[[1]]
    name = paste(name,collapse='__')
    #outpdf_name = paste(c(name,x[2]),collapse='__')
    #outpdf_name = paste(c(out,outpdf_name,'.png'),collapse='')
    #print(outpdf_name)
    #png(outpdf_name,width=1000,height=600,units='px',res=150)
    #png(outpdf_name)
        p = ggplot(subset(gene, molecule == x[1] & gene == x[2]),
        aes(xmin = start, xmax = end, y = strand)) +
        geom_gene_arrow(arrow_body_height = grid::unit(6, "mm"),arrowhead_width = grid::unit(7, "mm"), arrowhead_height = grid::unit(7, "mm")) +
		#geom_gene_arrow()+
        #geom_gene_label(aes(label = gene)) +
        geom_subgene_arrow(
        data = subset(sub_gene, molecule == x[1] & gene == x[2]),
        aes(xsubmin = from, xsubmax = to, fill = subgene),
		arrow_body_height = grid::unit(6, "mm")
        ) +
        geom_subgene_label(
        data = subset(sub_gene, molecule == x[1] & gene == x[2]),
        aes(xsubmin = from, xsubmax = to, label = subgene),
        min.size = 0
        )+
		labs(title=x[2])
		plot_continue(p,out)
        #print(plotpng)
    #dev.off()

}



enrich_last_plot <- function(term,feature){

	if(length(feature$y_method)!=0){
		newterm = term@result
		#newterm = newterm[-order(newterm[[feature$y_method]]),]
		#newterm = factor(newterm$Description,levels=unique(newterm$Description))
		if(length(feature$group)==0){
			
print(str(newterm))
			newterm = tail(newterm[order(-newterm[[feature$y_method]]),],10)
			print(newterm)
			newterm$Description = factor(newterm$Description,levels=unique(newterm$Description))
			print(newterm)
			p <- ggplot(newterm,aes(Description,-log10(get(feature$y_method)))) + geom_bar(stat='identity',,fill='red',position=position_dodge(2),width=0.4)  + labs(title=feature$plot_type,y=quote(-log[10](p-value)))+ coord_flip() 
		}
		else{
			print(str(newterm))
			newterm = tail(newterm[order(-newterm[[feature$y_method]],newterm[['ONTOLOGY']]),],10)
			newterm$Description = factor(newterm$Description,levels=unique(newterm$Description))
			p <- ggplot(newterm,aes(Description,-log10(get(feature$y_method)),fill=get(feature$group))) + geom_bar(stat='identity',width=0.4) + scale_fill_brewer(palette='Set1') + labs(title=feature$plot_type,y=quote(-log[10](p-value)),fill='Type') + coord_flip() 
		}
	}
	else{


	if(feature$plot_type=='enanalyse_go'){
		print(str(term))
		p <- barplot(term,showCategory=10,title='GO')
	}
	else{
		p <- barplot(term,showCategory=10,title='KEGG')
	}	

	}
	feature$mutiformat='bar'
	plot_continue(p,feature)

	p <- dotplot(term)
	feature$mutiformat='dot'
	plot_continue(p,feature)

	p <- cnetplot(term)
	feature$mutiformat='cnet'
	plot_continue(p,feature)

}

go_noref <- function(data,feature){
	table_path = paste0(feature$out,'/sample')
	table <- read.table(table_path,header = TRUE,sep = "\t",stringsAsFactors=F)
	print(str(table))
	
	genes <- data[[feature$gene_name]]
	gene_GOID <- data.frame(table[[feature$background_gkid]],table[[feature$background_gene]],stringsAsFactors=F)
	GO_names <- go2term(table[[feature$background_gkid]])
	enrich_go <- enricher(gene=genes,pvalueCutoff=feature$pvalue,pAdjustMethod=feature$method,minGSSize=10,maxGSSize=500,TERM2GENE=gene_GOID,TERM2NAME=GO_names)


	enrich_file = paste0(feature$out,'/go_enrich_noref_result.txt')
    write.table(enrich_go,enrich_file,sep='\t')

	enrich_last_plot(enrich_go,feature)
}	

go_ref <- function(data,feature){
	#library('org.Dr.eg.db')
	genes <- data[[feature$gene_name]]
	
	if(feature$gene_name_type == 'symbol'){
		enrich.go <- enrichGO(gene = genes,OrgDb = feature$ref_db,keyType='SYMBOL',ont = 'ALL',pAdjustMethod = feature$method,pvalueCutoff=feature$pvalue,qvalueCutoff = 0.2)
	}
	else if(feature$gene_name_type == 'ebi'){
		enrich.go <- enrichGO(gene = genes,OrgDb = feature$ref_db,keyType='ENSEMBL',ont = 'ALL',pAdjustMethod = feature$method,pvalueCutoff=feature$pvalue,qvalueCutoff = 0.2)
	}	
	else{
		enrich.go <- enrichGO(gene = genes,OrgDb = feature$ref_db,keyType='ENTREZID',ont = 'ALL',pAdjustMethod = feature$method,pvalueCutoff=feature$pvalue,qvalueCutoff = 0.2)
	}

	enrich_file = paste0(feature$out,'/go_enrich_ref_result.txt')
    write.table(enrich.go,enrich_file,sep='\t')
	

	enrich_last_plot(enrich.go,feature)
}

kegg_ref <- function(data,feature){
	library(AnnotationHub)		

	if(feature$ref_db[1]=='org.Dr.eg.db'){
		library(org.Dr.eg.db)
	}

	if(feature$ref_db[1]=='org.Hs.eg.db'){
		library(org.Hs.eg.db)
	}

	if(feature$ref_db[1]=='org.Mm.eg.db'){
		library(org.Mm.eg.db)
	}

	if(feature$ref_db[1]=='org.Rn.eg.db'){
		library(org.Rn.eg.db)
	}

	if(feature$ref_db[1]=='org.Dm.eg.db'){
		library(org.Dm.eg.db)
	}

	if(feature$ref_db[1]=='org.Ce.eg.db'){
		library(org.Ce.eg.db)
	}



if(feature$ref_db[1]=='org.At.tair.db'){
		library(org.At.tair.db)
	}

	if(feature$ref_db[1]=='org.Sc.sgd.db'){
		library(org.Sc.sgd.db)
	}

	if(feature$ref_db[1]=='org.Gg.eg.db'){
		library(org.Gg.eg.db)
	}

	if(feature$ref_db[1]=='org.Pt.eg.db'){
		library(org.Pt.eg.db)
	}

	if(feature$ref_db[1]=='org.Ss.eg.db'){
		library(org.Ss.eg.db)
	}

	genes <- data[[feature$gene_name]]

	if(feature$gene_name_type == 'symbol'){
		entrez_id = na.omit(mapIds(x=get(feature$ref_db[1]),keys=genes,keytype='SYMBOL',column='ENTREZID'))
	}
	else if(feature$gene_name_type == 'ebi'){
		entrez_id = na.omit(mapIds(x=get(feature$ref_db[1]),keys=genes,keytype='ENSEMBL',column='ENTREZID'))
	}
	else{
		entrez_id = na.omit(genes)
	}




	kegg <- enrichKEGG(gene=entrez_id,keyType='kegg',organism=feature$ref_db[2],pAdjustMethod=feature$method,pvalueCutoff=feature$pvalue,qvalueCutoff=0.2)

	enrich_file = paste0(feature$out,'/kegg_enrich_ref_result.txt')
    write.table(kegg,enrich_file,sep='\t')

	enrich_last_plot(kegg,feature)
}

kegg_noref <- function(data,feature){
	

	table_path = paste0(feature$out,'/sample')
	table <- read.table(table_path,header = TRUE,sep = "\t",stringsAsFactors=F)
	print(str(table))
	
	genes <- data[[feature$gene_name]]
	gene_KeID <- data.frame(table[[feature$background_gkid]],table[[feature$background_gene]],stringsAsFactors=F)
	ko_id <- bitr_kegg(table[[feature$background_gkid]], "kegg", "Path", "ko")
	ko_name <- ko2name(ko_id$Path)
	ko_k_name <- merge(ko_id,ko_name,by.ko_id ='Path',by.ko_name ='ko')
	k_name <- ko_k_name[c(1,4)]
	

	
	enrich_ke <- enricher(gene=genes,pvalueCutoff=feature$pvalue,pAdjustMethod=feature$method,minGSSize=10,maxGSSize=500,TERM2GENE=gene_KeID,TERM2NAME=k_name)

	enrich_file = paste0(feature$out,'/kegg_enrich_ref_result.txt')
    write.table(enrich_ke,enrich_file,sep='\t')

	enrich_last_plot(enrich_ke,feature)
}


enrich_analyse <- function(data,feature){
	library(clusterProfiler)
	if(feature$plot_type == 'enanalyse_go'){
		
		if(length(feature$ref_db)==0){
			
				go_noref(data,feature)
		}
		else{
			go_ref(data,feature)
		}
	}
	else{
		if(length(feature$ref_db)==0){
			kegg_noref(data,feature)		
		}
		else{
			kegg_ref(data,feature)
			
		}
	}
}

ppi_trans_id <- function(data,feature){
	
	library(AnnotationHub)		

	if(feature$ref_db[1]=='org.Dr.eg.db'){
		library(org.Dr.eg.db)
	}

	if(feature$ref_db[1]=='org.Hs.eg.db'){
		library(org.Hs.eg.db)
	}

	if(feature$ref_db[1]=='org.Mm.eg.db'){
		library(org.Mm.eg.db)
	}

	if(feature$ref_db[1]=='org.Rn.eg.db'){
		library(org.Rn.eg.db)
	}

	if(feature$ref_db[1]=='org.Dm.eg.db'){
		library(org.Dm.eg.db)
	}

	if(feature$ref_db[1]=='org.Ce.eg.db'){
		library(org.Ce.eg.db)
	}



	if(feature$ref_db[1]=='org.At.tair.db'){
		library(org.At.tair.db)
	}

	if(feature$ref_db[1]=='org.Sc.sgd.db'){
		library(org.Sc.sgd.db)
	}

	if(feature$ref_db[1]=='org.Gg.eg.db'){
		library(org.Gg.eg.db)
	}

	if(feature$ref_db[1]=='org.Pt.eg.db'){
		library(org.Pt.eg.db)
	}

	if(feature$ref_db[1]=='org.Ss.eg.db'){
		library(org.Ss.eg.db)
	}

	genes <- data[[feature$gene_name]]

	if(feature$gene_name_type == 'symbol'){
		entre_id = na.omit(select(x=get(feature$ref_db[1]),keys=genes,keytype='SYMBOL',column='ENTREZID'))
	}
	else if(feature$gene_name_type == 'ebi'){
		entre_id = na.omit(select(x=get(feature$ref_db[1]),keys=genes,keytype='ENSEMBL',column='ENTREZID'))
	}
	else{
		entre_id = na.omit(select(x=get(feature$ref_db[1]),keys=genes,keytype='ENTREZID',column='SYMBOL'))
	}
	entre_id = as.vector(entre_id)
	print(str(entre_id))
	string_db <- STRINGdb$new(version="10", species=as.numeric(feature$stringdb),score_threshold=feature$sth[1], input_directory="")
	data_mapped <- entre_id %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", removeUnmappedRows = TRUE)
	data_mapped_file = paste0(feature$out,'/idmapping')

	python_path_pre = strsplit(feature$out,"/")[[1]]
	print(python_path_pre)
	print(python_path_pre[1:(length(python_path_pre)-2)])
	python_path = c(python_path_pre[1:(length(python_path_pre)-2)],'ppi.py')
	
	python_sc = paste(python_path,collapse="/")
	print(python_sc)
	backgroud_file = paste(paste(python_path_pre[1:(length(python_path_pre)-2)],collapse="/"),"/stringdb/v11.5/",feature$stringdb,".*",sep='')
	comond_line = paste(python_sc,data_mapped_file,backgroud_file)
	print(comond_line)

    write.table(data_mapped[c('STRING_id','SYMBOL')],data_mapped_file,sep='\t',row.names=F,col.names=F,quote=FALSE)

	system(comond_line)

	out_file_name = data_mapped_file = paste0(feature$out,'/ppi_result.txt')
	ppi_data = read.delim(out_file_name,stringsAsFactors=F,header=T)
	print(str(ppi_data))
	return(ppi_data)
	#return(ensembl_id)
	
}

ppi_ggraph <- function(data,feature){
	#针对A->B,B->A的情况去除
	data_new <- t(apply(data,1,sort,decreasing=T))
	colnames(data_new) = colnames(data)
	data_new <- as.data.frame(data_new,stringsAsFactors=F)
	data_new$combined_score <- as.numeric(data_new$combined_score)
	data <- data_new %>% distinct()
	#针对A->B,B->A的情况去除

	data <- subset(data,combined_score>= feature$sth[1],drop=F)
	if(feature$removed == 'yes'){
		
		
		data <- data %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
  			mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
  			filter(!(from_c == 1 & to_c == 1)) %>%
  			dplyr::select(1,2,3)
	}
	
	print(data)
	print(str(data))
	data_note <- data %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
	net_ppi <- igraph::graph_from_data_frame(d=data,vertices=data_note,directed = F)
	igraph::V(net_ppi)$deg <- igraph::degree(net_ppi)
	igraph::V(net_ppi)$size <- igraph::degree(net_ppi)/5
	igraph::E(net_ppi)$width <- igraph::E(net_ppi)$combined_score/10
	if( length(feature$repel)!=0 && feature$repel=='yes'){
		feature$repel=T
	}else{
		feature$repel=F
	}
	if(feature$node_text=='yes'){
		if(feature$layout == 'linear'){
			p <- ggraph(net_ppi,layout = feature$layout, circular = TRUE)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+
  				geom_node_text(aes(filter=deg>feature$node_text_greater[1], label=name), size = feature$node_text_size[1], repel = feature$repel)+
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
		}else{
			if(feature$layout == 'centrality'){
			p <- ggraph(net_ppi,layout = feature$layout, cent = deg)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+
  				geom_node_text(aes(filter=deg>feature$node_text_greater[1], label=name), size = feature$node_text_size[1], repel = feature$repel)+
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
			}else{
				p <- ggraph(net_ppi,layout = feature$layout)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+
  				geom_node_text(aes(filter=deg>feature$node_text_greater[1], label=name), size = feature$node_text_size[1], repel = feature$repel)+
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
			}	  
		}
	}else{

		if(feature$layout == 'linear'){
			p <- ggraph(net_ppi,layout = feature$layout, circular = TRUE)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+  				
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
		}else{
			if(feature$layout == 'centrality'){
			p <- ggraph(net_ppi,layout = feature$layout, cent = deg)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+  				
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
			}else{
				p <- ggraph(net_ppi,layout = feature$layout)+
  				get(feature$edge_style)(aes(edge_width=width), color = feature$edge_color, show.legend = F)+
  				geom_node_point(aes(size=size), color=feature$node_color, alpha=feature$node_alpha[1])+  				
  				scale_edge_width(range = c(0.2,1))+
  				scale_size_continuous(range = c(1,10) )+
  				guides(size=F)+
  				theme_graph()
			}	  
		}

	}
	return(p)
}

ppi_analyse <- function(data,feature){
	library(tidyverse)
	library(clusterProfiler)
	library(STRINGdb)
	library(igraph)
	library(ggraph)
	ppi_data <- ppi_trans_id(data,feature)
	p <- ppi_ggraph(ppi_data,feature)
	plot_continue(p,feature)

}

deanalyse_pca <- function(dds,feature){
	

	rld <- rlog(dds)
	p <- plotPCA(rld, intgroup=c("condition"))+ stat_ellipse(type='t',level = 0.95,show.legend=F)+geom_text(aes(label=name),family='sans',size=3,vjust='top',check_overlap=T,hjust='center',nudge_y=-1,nudge_x=0,show.legend=F)+theme_light()
	feature$mutiformat='pca'
	plot_continue(p,feature)
}

deanalyse_pca_batch <- function(dds,feature){
	rld <- rlog(dds)
	p <- plotPCA(rld, intgroup=c("batch")) + stat_ellipse(type=feature$typed,level = as.numeric(feature$level),show.legend=F)
	data_new = plotPCA(rld, intgroup=c("batch"),returnData=T)
	print(str(data_new))
	plot_continue(p,feature)
}

normalize <- function(x){
  return ((x-mean(x))/(max(x)-min(x)))
}

deanalyse_valcannoone <- function(resdata,feature){
	library(ggrepel)
	print('valcannoone')
	resdata$selectedgene <- ifelse(resdata$change=='Down' | resdata$change=='Up',resdata[[1]],NA)
	print(str(resdata))
	valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
    geom_point(alpha=0.8, size=1) +
    theme_bw(base_size=15) +
    theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
    ) +
    ggtitle("") +
    #scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "Equally")) +
    geom_vline(xintercept=c(-abs(feature$fc), abs(feature$fc)), lty=2, col="gray", lwd=0.5) +
    #geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)
	geom_hline(yintercept=-log10(abs(feature$p)), lty=2, col="gray", lwd=0.5)+
	geom_text_repel(aes(label=selectedgene), color="black",size=2.5,box.padding=0.2,point.padding=0.3,segment.colour = "black",segment.size = 0.3,force = 1,show.legend = F)+theme_light()+
	scale_colour_brewer(palette='Accent')
	#ggrepel::geom_text_repel(aes(label=selectedgene),color="black")
	feature$mutiformat='valcano'
	plot_continue(valcano,feature)
}

denanlyse_radar <- function(resdata,table,feature){
	library(ggradar)
	end <- 8+nrow(table)-1
	if(nrow(resdata) < abs(feature$top_fc)){
		feature$top_fc = nrow(resdata)
	}
	resdata = na.omit(resdata)
	resdata <- resdata[order(abs(resdata$log2FoldChange),decreasing=TRUE,na.last=NA),]


	all_diff_gene_deseq2_expr_50 <- resdata[1:feature$top_fc,1:7]
	all_diff_gene_deseq2_expr_50[,2:7] = as.data.frame(apply(all_diff_gene_deseq2_expr_50[,2:7],2,normalize))
	

	print(str(all_diff_gene_deseq2_expr_50))
	p <- ggradar(all_diff_gene_deseq2_expr_50,grid.min = -1,
        grid.mid = 0, grid.max = 1,
		values.radar = c("-1", "0", "1"),
        gridline.min.colour = "grey",
        gridline.mid.colour = "blue", gridline.max.colour = "orange",
        axis.label.size = 3, axis.line.colour = "grey",
        legend.text.size = 8, legend.position = "bottom",
        background.circle.colour = "white",
        background.circle.transparency = 0.1,
		group.point.size = 3,
		group.line.width = 1,
        )
	p = p + theme_light() + theme(legend.position = "bottom")
	#p = p + theme_light(legend.position = "bottom")
	feature$mutiformat='radar'
	plot_continue(p,feature)


}

denanlyse_heatmap <- function(resdata,table,feature){
	library(pheatmap)
	end <- 8+nrow(table)-1

	if(nrow(resdata) < abs(feature$top)){
		feature$top = nrow(resdata)
	}

	all_diff_gene_deseq2_expr_50 <- resdata[1:feature$top,8:end]
	
	rownames(all_diff_gene_deseq2_expr_50) <- resdata[1:feature$top,]$Row.names
	print(str(all_diff_gene_deseq2_expr_50))
	#write.table(all_diff_gene_deseq2_expr_50,"/home/hli/Djangotest/GBTP/static/magicRversion1/data/heatmap",sep='\t',row.names=T)
	print(str(table))
	print('dasdasdas')
	print(table)
	#p<- pheatmap(all_diff_gene_deseq2_expr_50, cluster_rows=TRUE, scale="row", annotation_col=table)
	feature$mutiformat='heat'
	feature$r_method = 'ggplot2'
	
	heatmap_all(all_diff_gene_deseq2_expr_50,feature,table)
	#plot_continue(p,feature)
}

denanlyse_batch <- function(data,table,feature){
	library("DESeq2")
	library("BiocParallel")
	#register(MulticoreParam(4))
	condition <- table[,2]
	batch <- table[,1]
	colnames(table) = c('batch','condition')
	print(table)
	print(condition)
	print(str(table))
	print(str(data))
	#构建样品信息矩阵

	multiple = max(data)/.Machine$integer.max
	if(multiple>1){
		data = round(data/ceiling(multiple))
		
	}

	if(feature['batch'] == 'DEseq'){
		print('batch');print(batch);print(condition)
		#table$condition = factor(table$condition,levels=unique(table$condition))
		#table$batch = factor(table$batch,levels=unique(table$batch))
		print(str(table))
		#factor(data[[feature$group]],levels=unique(data[[feature$group]]))
		data <- DESeqDataSetFromMatrix(data, colData=table, design= ~ condition+batch)
	}else{
		data <- DESeqDataSetFromMatrix(data, colData=table, design= ~ condition)
	}
	#构建dds矩阵
	dds <- data[ rowSums(counts(data)) > 1, ]
	#过滤count数都为0的数据
	
	
	deanalyse_pca_batch(dds,feature)
	
}

deanalyse <- function(data,table,feature){
	library("DESeq2")
	library("BiocParallel")
	register(MulticoreParam(16))
	condition <- table[,2]
	batch <- table[,1]
	colnames(table) = c('batch','condition')
	print(table)
	print(condition)
	print(str(table))
	print(str(data))
	#构建样品信息矩阵

	multiple = max(data)/.Machine$integer.max
	if(multiple>1){
		data = round(data/ceiling(multiple))
		
	}

	if(feature['batch'] == 'DEseq'){
		print('batch');print(batch);print(condition)
		data <- DESeqDataSetFromMatrix(data, colData=table, design= ~ condition+batch)
	}else{
		data <- DESeqDataSetFromMatrix(data, colData=table, design= ~ condition)
	}
	#构建dds矩阵
	
	check_length_f = paste0(feature['out'],'/mysql_selected_length.txt')
	if (file.exists(check_length_f)){
		length_s <- read.table(check_length_f, header=T, row.names=1, com='', quote='',check.names=F, sep="\t")
		print(str(length_s))                       
		data = DESeq2::estimateSizeFactors(data, normMatrix = length_s)  # 长度校准
	}
	print('Size end')
	

	keep <- rowSums(counts(data) >= feature$reads_number_de) >= feature$filter_samples_de
	dds <- data[keep, ]
	#dds <- data[ rowSums(counts(data)) > 1, ]
	#过滤count数都为0的数据
	mean_dds = dds
	if (length(intersect(c('deanalyse_pca','deanalyse_all'),feature$format))!=0){ 
		deanalyse_pca(dds,feature)
	}
	tryCatch({
		dds <- DESeq(dds,parallel = TRUE)

		resultsNames(dds)
	res <- results(dds)
	#获取结果并赋值给res
	res <- res[order(res$padj),]
	# 按照padj的大小将res重新排列
	resdata <- merge (as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
	#画火山图	
	
	print(feature)

	resdata$change <- as.factor(
    	ifelse(
        	resdata$padj<=(feature$p) & abs(resdata$log2FoldChange)>=abs(feature$fc),
        	ifelse(resdata$log2FoldChange>=abs(feature$fc), "Up", "Down"),
        	"Equally"
    	)
	)

	all_file = paste0(feature$out,'/all_gene_result.txt')
    write.table(resdata,all_file,sep='\t',row.names=F,quote=FALSE)
	
	all_file = paste0(feature$out,'/up_gene_result.txt')
    write.table(subset(resdata,change=='Up',drop=F),all_file,sep='\t',row.names=F,quote=FALSE)
	
	all_file = paste0(feature$out,'/down_gene_result.txt')
    write.table(subset(resdata,change=='Down',drop=F),all_file,sep='\t',row.names=F,quote=FALSE)

	if (length(intersect(c('deanalyse_valcannoone','deanalyse_all'),feature$format))!=0){ 
		deanalyse_valcannoone(resdata,feature)
	}

	if(feature$de_module_fc == 'up'){
		resdata_fc = subset(resdata,change=='Up',drop=F)
	}
	else if(feature$de_module_fc == 'down'){
		resdata_fc = subset(resdata,change=='Down',drop=F)
	}
	else{
		resdata_fc = resdata
		print('i am all_fc')
	}

	if (length(intersect(c('denanlyse_heatmap','deanalyse_all'),feature$format))!=0){ 
		

		denanlyse_radar(resdata_fc,table,feature)
	}

	if(feature$de_module == 'up'){
		resdata = subset(resdata,change=='Up',drop=F)
	}
	else if(feature$de_module == 'down'){
		resdata = subset(resdata,change=='Down',drop=F)
	}
	else{
		print('i am all')
	}



	if (length(intersect(c('denanlyse_heatmap','deanalyse_all'),feature$format))!=0){ 
		denanlyse_heatmap(resdata,table,feature)

		
	}

	



	},error=function(e){
		mean_dds <- DESeq(mean_dds,parallel = T,fitType="mean")
		
		print("DESeq2::DESeq(,fitType='mean'). The fitting method will only adopt 'mean' when parametric, local and local all fail, which means that your data seems to have too many low-expression genes.")
		resultsNames(mean_dds)
	res <- results(mean_dds)
	#获取结果并赋值给res
	res <- res[order(res$padj),]
	# 按照padj的大小将res重新排列
	resdata <- merge (as.data.frame(res),as.data.frame(counts(mean_dds,normalize=TRUE)),by="row.names",sort=FALSE)
	#画火山图	
	
	print(feature)

	resdata$change <- as.factor(
    	ifelse(
        	resdata$padj<=(feature$p) & abs(resdata$log2FoldChange)>=abs(feature$fc),
        	ifelse(resdata$log2FoldChange>=abs(feature$fc), "Up", "Down"),
        	"Equally"
    	)
	)

	all_file = paste0(feature$out,'/all_gene_result.txt')
    write.table(resdata,all_file,sep='\t',row.names=F,quote=FALSE)
	
	all_file = paste0(feature$out,'/up_gene_result.txt')
    write.table(subset(resdata,change=='Up',drop=F),all_file,sep='\t',row.names=F,quote=FALSE)
	
	all_file = paste0(feature$out,'/down_gene_result.txt')
    write.table(subset(resdata,change=='Down',drop=F),all_file,sep='\t',row.names=F,quote=FALSE)

	if (length(intersect(c('deanalyse_valcannoone','deanalyse_all'),feature$format))!=0){ 
		deanalyse_valcannoone(resdata,feature)
	}

	if(feature$de_module_fc == 'up'){
		resdata_fc = subset(resdata,change=='Up',drop=F)
	}
	else if(feature$de_module_fc == 'down'){
		resdata_fc = subset(resdata,change=='Down',drop=F)
	}
	else{
		resdata_fc = resdata
		print('i am all_fc')
	}

	if (length(intersect(c('denanlyse_heatmap','deanalyse_all'),feature$format))!=0){ 
		

		denanlyse_radar(resdata_fc,table,feature)
	}



	if(feature$de_module == 'up'){
		resdata = subset(resdata,change=='Up',drop=F)
	}
	else if(feature$de_module == 'down'){
		resdata = subset(resdata,change=='Down',drop=F)
	}
	else{
		print('i am all')
	}

	if (length(intersect(c('denanlyse_heatmap','deanalyse_all'),feature$format))!=0){ 
		denanlyse_heatmap(resdata,table,feature)
		
		
	}
	return(message("错误，数据是三列"))
	})
	


}





pheatmap_plot <- function(data,feature,group){
	library(pheatmap)
	library(ggplotify)
	if(feature$legend == 'yes'){
		feature$legend=T
	}else{
		feature$legend=F
	}

	if(feature$number == 'no'){
		feature$number = F
	}else{
		feature$number=T
	}

	if(feature$title == ''){
		feature$title=NA
	}

	if(feature$process != 'scale' ){
		feature$process='none'
	}
	else{
		feature$process='row'
	}
	
	if(length(feature$distance)!=0){
		if(length(feature$sample)!=0){
			if(length(feature$median) == 0){
				feature$median = 'complete'
			}
			
			
			
			tablefile = paste0(feature$out,'/sample')
			#print(tablefile)
			table = read.delim(tablefile,stringsAsFactors=F,row.name=1)
			#print(table)
			#print(colnames(data))
			#print(str(data))
			p <- pheatmap(data, cluster_rows=T,clustering_distance_rows=feature$distance,clustering_method=feature$method,annotation_col=table,scale=feature$process,legend=feature$legend,display_numbers=feature$number,main=feature$title,angle_col=feature$angle,fontsize=feature$size)
		}
		else{	
			p <- pheatmap(data, cluster_rows=T,clustering_distance_rows=feature$distance,clustering_method=feature$method,scale=feature$process,legend=feature$legend,display_numbers=feature$number,main=feature$title,angle_col=feature$angle,fontsize=feature$size)
		}
	}
	else{
		if(length(feature$sample)!=0){
			#print(feature$out)
			tablefile = paste0(feature$out,'/sample')
			#print(tablefile)
			table = read.delim(tablefile,stringsAsFactors=F,row.name=1)
			#print(str(table))
			p <- pheatmap(data, cluster_rows=T,annotation_col=table,scale=feature$process,legend=feature$legend,display_numbers=feature$number,main=feature$title,angle_col=feature$angle,fontsize=feature$size)

		}		
		else{	
			print('hres')
			p <- pheatmap(data, cluster_rows=F,cluster_cols=F,legend=feature$legend,display_numbers=feature$number,scale=feature$process,main=feature$title,angle_col=feature$angle,fontsize=feature$size)
		}	
	}
	p <- as.ggplot(p)
	return (p)

}

heatmap_pheatmap <- function(data,feature){
	if(feature$process !='0' ){
		data <- heatmap_process(data,feature)
	}	
	group = 0
	

	p <- pheatmap_plot(data,feature,group)	
	plot_continue(p,feature)

}

log10del <- function(x){
	x <- x+1
	return(log10(x))
}

heatmap_process <- function(df,feature){
	library(reshape2)
	library(ggtree)
	library(scales)
		if(feature$process == 'log10'){
			df <- apply(df,1,log10del)
			df <- as.data.frame(t(df),stringsAsFactors=FALSE)
			return(df)
		}
		if(feature$process == 'scale' && feature[['r_method']] == 'pheatmap'){
			return (df)
		}

		df <- apply(df,1,get(feature$process))
		df <- as.data.frame(t(df),stringsAsFactors=FALSE)
		return (df)
	
}

counstruct_group <- function(feature,table=0){
	#library(dplyr)
	#samplefile <- paste0(feature$out,'/sample')
	#samplefile <- paste0(feature$out,'/mysql_condition.txt')
	#sample_data <- read.delim(samplefile,stringsAsFactors=F)
	sample_data = table
	groups <- unique(sample_data[[2]])
	rowlen = length(rownames(sample_data))
	sample_data$p = replicate(rowlen,colnames(sample_data)[2])
	sample_data$batch = rownames(sample_data)
	temp_second_col <- colnames(sample_data)[2]
	print(temp_second_col)
	colnames(sample_data) = c('sample','condition','p')
	print('asdasd')
	print(str(sample_data))

	group <- ggplot(sample_data,aes(sample,y=p,fill=condition)) + geom_tile(color='black') + scale_y_discrete(position="right",labels=c(temp_second_col),expand=c(0,0)) + theme_minimal() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + labs(fill=temp_second_col)
	
	return (group)

}

counstruct_tree <- function(data,feature){
	library(ggtree)
	
	gg_row <- ggtree(hclust(dist(data,method=feature$distance),method=feature$method),layout="rectangular",branch.length="none")
	gg_col <-  ggtree(hclust(dist(t(data),method=feature$distance),method=feature$method)) + layout_dendrogram()
	return_p = list(a=gg_col,b=gg_row)
	return(return_p)
}

heatmap_plot <- function(df,feature,group,gg_tree){

	library(reshape2)
	library(ggtree)
	library(scales)
	library(aplot)
	#library(ggplotify)
	#library(patchwork)
	
	print(str(df))
	rowlen = length(rownames(df))
	collen = length(colnames(df))
	df <- cbind(name=row.names(df),df,stringsAsFactors=FALSE)

	df <- melt(df,id.vars='name')
	df$variable = as.vector(df$variable)
	print(str(df))
	#df$x <- rep(c(1:collen),each=rowlen)
	#df$y <- rep(c(1:rowlen),collen)		
	print(str(df))
	
	#heatmap_p <- ggplot(df,aes(variable,name,fill=value)) + geom_raster() + scale_fill_gradient2(low="#003366", high="#990033", mid="white") + geom_tile(color='black') + theme_minimal() + xlab(NULL) + ylab(NULL)
	heatmap_p <- ggplot(df,aes(variable,name))+scale_color_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))+theme_bw()+geom_point(aes(size=value,color=value))+  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =45,hjust =1))+xlab(NULL) + ylab(NULL)+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+geom_vline(xintercept=rep(seq(1.5, 1.5+collen-2), 1),size=.2)#+geom_vline(xintercept=c(3.5),size=.5)


	

	if(length(feature$distance)!=0){
		if(length(feature$sample)!=0){
			samplefile <- paste0(feature$out,'/sample')
	#samplefile = paste0("/home/hli/Djangotest/GBTP/static/magicRversion1/data/",feature$out,"/sample")
			#sample_data <- read.delim(samplefile,stringsAsFactors=F)
			if(feature$hclust == 'both'){
				p <- heatmap_p %>% insert_top(group, height = .05) %>% insert_top(gg_tree$a,height=.1)  %>% insert_left(gg_tree$b,width=.2)
			}
			if(feature$hclust == 'row'){
				p <- heatmap_p %>% insert_top(group, height = .05) %>% insert_left(gg_tree$b,width=.2)
			}
			if(feature$hclust == 'col'){
				p <- heatmap_p %>% insert_top(group, height = .05) %>% insert_top(gg_tree$a,height=.1)
			}
		}
		else{
			if(feature$hclust == 'both'){
				p <- heatmap_p %>% insert_left(gg_tree$b,width=.2) %>% insert_top(gg_tree$a,height=.2)
			}
			if(feature$hclust == 'row'){
				p <- heatmap_p %>% insert_left(gg_tree$b,width=.2)
			}
			if(feature$hclust == 'col'){
				p <- heatmap_p %>% insert_top(gg_tree$a,height=.2)
			}

			#p <- heatmap_p %>% insert_left(gg_tree$b,width=.2) %>% insert_top(gg_tree$a,height=.2)
			#p <- as.ggplot(p)
			#print(str(p))
			#p <- gg_tree$b + heatmap_p
		}
	}
	else{
		
		if(length(feature$sample)!=0){
			#p <- group
			p <- heatmap_p %>% insert_top(group, height = .05)
			#p <- as.ggplot(p)
			#p <- group / plot_spacer() /heatmap_p + plot_layout(heights=c(1,-1.6,35),guides = 'collect')
		}else{
			p <- heatmap_p
		}			
		
		#scale_fill_gradient(low='steelblue',high='red')
		#scale_fill_gradient2(low=muted('blue'),midpoint=0.5,high=muted('red'),mid='yellow')
		#scale_fill_gradient(low='blue',high='red')
		#scale_fill_gradient2(low=muted('blue'),midpoint=(max(df$value)+min(df$value))/2,high=muted('red'),mid='white')
		#scale_fill_gradient2(low=muted('green'),high=muted('red'),mid='white',midpoint=(max(df$value)+min(df$value))/2)
		#scale_fill_gradient2(low=muted('blue'),midpoint=median(df$value),high=muted('red'),mid='white')
		#scale_fill_gradientn(values = seq(-2,2,0.2),colours = c('cyan','blue','green','orange','red'))
		#scale_fill_gradient2(low=muted('blue'),midpoint=median(df$value),high=muted('red'))
		#scale_fill_gradient(low='blue',high='red')
		#scale_fill_gradient2(low=muted('blue'),midpoint=0,high=muted('red'))		
	}
	return (p)
}

heatmap_all <- function(data,feature,table=0){
	if(feature$process !='0'){
		data_old <- heatmap_process(data,feature)
		colnames(data_old) = colnames(data)
		data = data_old
	}	
	group = 0
	if(length(feature$sample)!=0){
		
		group <- counstruct_group(feature,table)
	}
	gg_tree = 0
	if(length(feature$distance)!=0){
		gg_tree <- counstruct_tree(data,feature)
	}

	p <- heatmap_plot(data,feature,group,gg_tree)	
	#p <- group
	plot_continue(p,feature)	
}


user_css <- function(p,conf){
	if(length(conf[[1]]$theme)!=0){
		if(conf[[1]]$theme %in% c("theme_tufte","theme_solarized","theme_excel","theme_hc")){
			library(ggthemes)			
		}		
		p <- p + get(conf[[1]]$theme)()
		conf = conf[-1]
	}
	#p <- p + get(conf[[1]]$theme)()
	#ggsave("scale.png",plot=p)
	#conf = conf[-1]
	for (css in conf){
		css <- construct_theme(css)	
		#print(css)	
		p <- p+css
	}
	return(p)
}

counstruct_text <- function(item){
#b <- list(plot.title=list(colour="red",hjust=NULL,vjust=NULL,size=NULL,margin=NULL))
	temp= item[[1]]
	a <- list(family=temp$family,face=temp$face,colour=temp$colour,size=temp$size,hjust=temp$hjust,vjust=temp$vjust,angle=temp$angle,lineheight=temp$lineheight,margin=temp$margin,debug=temp$debug)
	item[[1]] = a
	return(item)
}

construct_line <- function(item){
	temp = item[[1]]
	a <- list(colour=temp$colour,size=temp$size,linetype=temp$linetype)
	item[[1]] = a
	return(item)
}

construct_rect <- function(item){
    temp = item[[1]]
    a <- list(colour=temp$colour,size=temp$size,fill=temp$fill)
    item[[1]] = a
    return(item)
}



counstruct_margin <- function(item){
	temp_int = item[[1]]$rbl
	temp_unit = item[[1]]$unit
	a <- get(item$subclass[1])(temp_int[1],temp_int[2],temp_int[3],temp_int[4],unit=temp_unit)
	item[[1]]=a
	item$subclass=NULL
	attr(item,"class")=c("theme","gg")
	return(item)
}

construct_ticks <- function(item){
	temp_int = item[[1]]$number
	temp_unit = item[[1]]$unit
	a <- unit(temp_int,temp_unit)
	item[[1]] = a
	item$subclass=NULL
	attr(item,"class")=c("theme","gg")
	return(item)	
}

construct_theme <- function(item){
	#if(names(item)[1] %in% c('axis.ticks.length.x','axis.ticks.length.y')){
    #    item <- construct_ticks(item)
    #    return(item)
    #}

	if(length(names(item))==1){
        attr(item,"class")=c("theme","gg")
        return(item)
    }

    if(item$subclass[1] == "margin"){
        item <- counstruct_margin(item)
        return(item)
    }

	if(item$subclass[1] == "unit"){
        item <- construct_ticks(item)
        return(item)
    }

	if(item$subclass[[2]][1] == "element_text" ){
		item <- counstruct_text(item)		
	}

	if(item$subclass[[2]][1] == "element_line"){
		item <- construct_line(item)
	}

	if(item$subclass[[2]][1] == "element_rect"){
        item <- construct_rect(item)
    }	

	attr(item,"class")=c("theme","gg")
	#print(item[[names(item)[1]]])
	attr(item[[names(item)[1]]],item$subclass[[1]]) = item$subclass[[2]]
	item$subclass = NULL
	print(str(item))
	return(item)
}

construct_continuous <- function(cp){
	subcp <- cp[[1]]
	if(length(subcp$breaks)==0){
		subcp$breaks = waiver()		
	}
	if (length(subcp$labels)==0){
		subcp$labels = waiver()
	}	
	if(length(subcp$expand)==0){
		subcp$expand = waiver()
	}
	if (length(subcp$limits)==0){
		subcp$limits = NULL
	}
	if (length(subcp$position)==0){
		if(names(cp) == "scale_x_continuous" || names(cp) == "scale_x_discrete"){
			subcp$position = "bottom"
		}
		else{
			subcp$position = "left"
		}
	}
	if (length(subcp$trans)==0){
		subcp$trans = "identity"
	}
	return(subcp)
}

user_scale <- function(p,conf){
	for(cp in conf){
		if (names(cp) %in% c("scale_x_continuous","scale_y_continuous")){
			subcp <- construct_continuous(cp)
			p <- p+get(names(cp))(breaks=subcp$breaks,labels=subcp$labels,position=subcp$position,expand = subcp$expand,limits = subcp$limits,trans = subcp$trans)
			next
		}
		
		if (names(cp) %in% c("scale_x_discrete","scale_y_discrete")){
			subcp <- construct_continuous(cp)
			p <- p+get(names(cp))(breaks=subcp$breaks,labels=subcp$labels,position=subcp$position,expand = subcp$expand)
			next
		}

		if (names(cp) == "labs"){
			temp <- cp$labs
			attr(temp,'class')="labels"
			print(temp)
			p <- p+temp	
			next
		}

        if(names(cp) == "coord_cartesian"){
            #if (length(cp$coord_cartesian$xlim) == 0){
			#	cp$coord_cartesian$xlim = NULL
            #}
            #if (length(cp$coord_cartesian$ylim) == 0){
			#	cp$coord_cartesian$ylim = NULL
            #}
			p <- p+coord_cartesian(xlim=cp$coord_cartesian$xlim,ylim=cp$coord_cartesian$ylim)
            next
        }

		if(names(cp) == "geom_text"){
			if(length(cp$geom_text$vjust) == 0){
				cp$geom_text$vjust="middle"
			}
			if(length(cp$geom_text$hjust) == 0){
				cp$geom_text$hjust="center"
			}
			if(length(cp$geom_text$nudge_y) == 0){
				cp$geom_text$nudge_y=0
			}
			if(length(cp$geom_text$nudge_x) == 0){
				cp$geom_text$nudge_x=0
			}
			if(length(cp$geom_text$size) == 0){
				cp$geom_text$size=4
			}
			if(length(cp$geom_text$family) == 0){
				cp$geom_text$family="sans"
			}
			if(length(cp$geom_text$check_overlap) == 0){
				cp$geom_text$check_overlap = F
			}
			else{
				cp$geom_text$check_overlap = T
			}

			temp=cp[[1]]
			#NA可以消除text
			#print(temp)
			p <- p+geom_text(aes(label=text),family=temp$family,size=temp$size,vjust=temp$vjust,check_overlap=temp$check_overlap,hjust=temp$hjust,nudge_y=temp$nudge_y,nudge_x=temp$nudge_x,show.legend=F)	
			next
		}

		

		if(names(cp) == "geom_dl"){
			if(length(cp$geom_dl$label)==0){
				cp$geom_dl$label="group"
			}
			if(length(cp$geom_dl$method)==0){
				cp$geom_dl$method="smart.grid"
			}
			library(directlabels)
			print(cp$geom_dl$label)
			a <- cp$geom_dl$label		
			p <- p+directlabels::geom_dl(aes(label=get(a)),method=cp$geom_dl$method)
			next
		}

		if(names(cp) == "annotate"){
			#if(length(cp$geom_dl$label)==0){
			#	cp$geom_dl$label=""
			#}
			if(length(cp$annotate$size)==0){
				cp$annotate$size=4
			}
			if(length(cp$annotate$colour)==0){
				cp$annotate$colour="black"
				
			}

			p <- p+ annotate('text',x=cp$annotate$x,y=cp$annotate$y,label=cp$annotate$label,size=cp$annotate$size,colour=cp$annotate$colour)			
			next
		}

		if(names(cp) %in% c("scale_fill_brewer","scale_colour_brewer")){
			#print(cp[[1]]$palette)
			p <- p + get(names(cp))(palette=(cp[[1]]$palette))
			#save(p,file='/home/hli/Djangotest/GBTP/static/magicRversion1/data/GEfSLmRuyO8d163gTeP/debug')
			next
		}
		
		
		if(names(cp) == 'scale_fill_gradient'){
			print(names(cp))
			p <- p + get(names(cp))(low=cp[[1]]$low,high=cp[[1]]$high) + guides(fill=guide_colorbar(reverse=T))
			#save(p,file='/home/hli/Djangotest/GBTP/static/magicRversion1/data/GEfSLmRuyO8d163gTeP/debug')
			next
		}

		if(names(cp) == 'scale_colour_gradient'){
			p <- p + get(names(cp))(low=cp[[1]]$low,high=cp[[1]]$high) + guides(color=guide_colorbar(reverse=T)) + scale_size()
			#save(p,file='/home/hli/Djangotest/GBTP/static/magicRversion1/data/GEfSLmRuyO8d163gTeP/debug')
			next
		}

		if(names(cp) == 'scale_size'){
			p <- p + scale_size(range=c(cp[[1]][1],cp[[1]][2]))
			#save(p,file='/home/hli/Djangotest/GBTP/static/magicRversion1/data/GEfSLmRuyO8d163gTeP/debug')
			next
		}



		if(names(cp) == "scale_fill_grey"){
			if(length(cp[[1]]$start) == 0){
				p <- p + scale_fill_grey()
			}
			else{
				p <- p + scale_fill_grey(start=cp[[1]]$start,end=cp[[1]]$end)
			}
			next
		}
		
		if(names(cp) %in% c("scale_fill_hue","scale_colour_hue")){
			if(length(cp[[1]]$h) == 0){
				cp[[1]]$h = c(0,360) + 15
			}
			if(length(cp[[1]]$c) == 0){
				cp[[1]]$c = 100	
			}
			if(length(cp[[1]]$l) == 0){
				cp[[1]]$l = 65
			}
			p <- p + get(names(cp))(h=cp[[1]]$h,c=cp[[1]]$c,l=cp[[1]]$l)
			next
		}

		if(names(cp) %in% c("scale_fill_manual","scale_colour_manual")){
			library(wesanderson)
			if (length(cp[[1]]$wes_palette) == 0){
				p <- p+get(names(cp))(values=cp[[1]]$values)		
				next
			}
			else{
				
				print(wes_palette('Zissou'))
				p <- p+get(names(cp))(values = wes_palette(cp[[1]]$wes_palette))
				next
			}
		}
		if(names(cp) %in% c('geom_point')){
			p <- p + get(names(cp))(size=cp[[1]]$size)
		}

		if(names(cp) %in% c('coord_fixed','coord_flip')){
			p <- p+get(names(cp))()
		}

		if(names(cp) %in% c('geom_signif')){
			library(ggsignif)
			if(cp[[1]]$asterisk == 'no'){
				cp[[1]]$asterisk=F
			}
			else{
				cp[[1]]$asterisk=T
			}
			compaired = combn(cp[[1]]$groups,2,simplify=F)
			print(cp[[1]])
			#
			p <- p + geom_signif(comparisons = compaired,map_signif_level=cp[[1]]$asterisk,test=cp[[1]]$method,size=cp[[1]]$size,textsize=cp[[1]]$textsize,y_position=cp[[1]]$position,annotations=cp[[1]]$annote)
		}


	}
	return(p)		
}

detail_plot <- function(plot.conf){
   	user.conf <- fromJSON(file=plot.conf)
	if(length(user.conf$flage_type)==0){
		plotitem_pathy <- paste0(user.conf$out,"/plot_",user.conf$count+user.conf$checkedbutton)
	}
	else{
		plotitem_pathy <- paste0(user.conf$out,"/plot_",user.conf$count+user.conf$checkedbutton,'_',user.conf$flage_type)
	}
   	
   	load(plotitem_pathy)#默认为p
	if((length(user.conf$plottype_temp)!=0 && user.conf$plottype_temp == 'heatmap' && attr(p,'class')[1]!='gg') || (length(user.conf$flage_type)!=0 && user.conf$flage_type == 'heat' && attr(p,'class')[1]!='gg')){
		library('aplot')
		print('hello')
		#p_ori <- p
		p_main <- p$plotlist[[1]]
		p_new <- user_scale(p_main,user.conf$scale)
		p_new_new <- user_css(p_new,user.conf$css)
		p$plotlist[[1]] <- p_new_new
		plot_continue(p,user.conf)
	}else{
   		p <- user_scale(p,user.conf$scale)
	#ggsave("scale.png",plot=p)
   		p <- user_css(p,user.conf$css)
	#print(str(p))
	#ggsave("scale.png",plot=p)
	#if((length(user.conf$plottype_temp)!=0 && user.conf$plottype_temp == 'heatmap' && attr(p,'class')[1]!='gg') || (length(user.conf$flage_type)!=0 && user.conf$flage_type == 'heat' && attr(p,'class')[1]!='gg')){
	#	print('hello02')
		
	#	p_ori$plotlist[[1]] <- p
	#	p <- p_ori
	#}

   		plot_continue(p,user.conf)
	}
}

plot_continue <- function(p,conf){
	if(length(conf$maxdetail)==0){
		if(length(conf$mutiformat) == 0){
        	main_file = paste0(conf$out,"/ordination_",conf$count+1,".png")#这个文件名称值得商榷
    #print(main_file)
			if (length(conf$sample_name)!=0){
				
				print(conf['batch'])
				if(conf['batch']=="DEseq"){
					p <- p+geom_text(aes(label=name),family='sans',size=3,vjust='top',check_overlap=T,hjust='center',nudge_y=-1,nudge_x=0,show.legend=F)
				}else{
					p <- p+geom_text(aes(label=text),family='sans',size=3,vjust='top',check_overlap=T,hjust='center',nudge_y=-1,nudge_x=0,show.legend=F)	
				}


			}
    		ggsave(main_file,plot=p)
    		save_plot <- paste0(conf$out,"/plot_",conf$count+1)
    		save(p,file=save_plot)

			save_plot <- paste0(conf$out,"/plot_",conf$count+1,".RData")
    		save(p,file=save_plot)
		}
		else{
			main_file = paste0(conf$out,"/ordination_",conf$count+1,'_',conf$mutiformat,".png")
			ggsave(main_file,plot=p)
			save_plot <- paste0(conf$out,"/plot_",conf$count+1,'_',conf$mutiformat)
			save(p,file=save_plot)

			save_plot <- paste0(conf$out,"/plot_",conf$count+1,'_',conf$mutiformat,".RData")
			save(p,file=save_plot)
		}	
    }
	else{
		#main_file = paste0(conf$out,"/ordination_",conf$count+conf$maxdetail,".png")#这个文件名称值得商榷	
		#ggsave(main_file,plot=p)
		#save_plot <- paste0(conf$out,"/plot_",conf$count+conf$maxdetail)
		#save(p,file=save_plot)
		print("i am here")
		if(length(conf$download$width)==0){
			conf$download$width=NA
		}
		if(length(conf$download$height)==0){
            conf$download$height=NA
        }
		if(length(conf$download$res)==0){
			conf$download$res=300
		}
		
		if(length(conf$flage_type)==0){
			main_file = paste0(conf$out,"/ordination_",conf$count+conf$maxdetail,".png")#这个文件名称值得商榷
			save_plot <- paste0(conf$out,"/plot_",conf$count+conf$maxdetail)
			save_plot_fc <- paste0(conf$out,"/plot_",conf$count+conf$maxdetail,".RData")
		}
		else{
			main_file = paste0(conf$out,"/ordination_",conf$count+conf$maxdetail,'_',conf$flage_type,".png")
			save_plot <- paste0(conf$out,"/plot_",conf$count+conf$maxdetail,'_',conf$flage_type)
			save_plot_fc <- paste0(conf$out,"/plot_",conf$count+conf$maxdetail,'_',conf$flage_type,".RData")
		}	

		
		
		
        ggsave(main_file,plot=p,width=conf$download$width,height=conf$download$height,dpi=conf$download$res)
		
       
        save(p,file=save_plot)
		save(p,file=save_plot_fc)


		for (type in conf$download$format){	
			if(length(conf$flage_type)==0){
				main_file = paste0(conf$out,"/ordination_",conf$count+conf$maxdetail,".",type)#这个文件名称值得商榷
		#print(main_file)	
			}else{
				main_file = paste0(conf$out,"/ordination_",conf$count+conf$maxdetail,'_',conf$flage_type,".",type)
			}	
			ggsave(main_file,width=conf$download$width,height=conf$download$height,dpi=conf$download$res,plot=p)
		}
	}
}

err_min <- function(x){
    mean(x)-sd(x)
}

err_max <- function(x){
    mean(x)+sd(x)
}


barmean_plot <- function(df,feature,comparison){
	
	library(ggsignif)
	nam_all <- names(df)
	name_leave = nam_all[which(nam_all!=feature$group)]
	p <- ggplot(data=df,aes(get(feature$group),get(name_leave[1])))
	if (feature$fill == '--- no mapping ---'){
		p <- p + stat_summary(geom="bar",fun=mean,color='black') 
	}
	else{
		p <- p + stat_summary(geom="bar",fun=mean,aes(fill=get(feature$group)),color='black')
		#p <- p + stat_summary(geom="bar",fun=mean,aes(linetype=get(feature$group))) 
	}

	if(feature$err_bar != 'no'){
		#p <- p + stat_summary(geom="errorbar",fun.min=err_min,fun.max=err_max,width=0.2)
		p <- p + stat_summary(geom="errorbar",fun.min=mean,fun.max=err_max,width=0.2)
	}


	if(feature$markers == 'no'){
		p <- p + labs(x=feature$group,y=name_leave,fill=feature$group)
		return (p)
	}
	else{
		#保持设计思路的完整性
		if(feature$asterisk == 'no'){
			feature$asterisk=F
		}
		else{
			feature$asterisk=T
		}
		#print(feature)
		#print(comparison)
		p <- p + geom_signif(comparisons=comparison,test=feature$smethod,step_increase=feature$interval,map_signif_level =feature$asterisk,size=feature$size,textsize=feature$textsize)
		p <- p + labs(x=feature$group,y=name_leave,fill=feature$group)
		return(p)		

	}

}

boxmean_plot <- function(df,feature,comparison){
	library(ggpubr)
	nam_all <- names(df)
	name_leave = nam_all[which(nam_all!=feature$group)]
	#palette aaas,lancet not yet test
	#best ucscgb
	palette_multicolor = "ucscgb"
	if (feature$fill == 'no' && feature$err_bar=='no'){
		p <- ggboxplot(df,x=(feature$group),y=(name_leave),palette=palette_multicolor) 
	}
	else if(feature$fill != 'no' && feature$err_bar=='no'){
		p <- ggboxplot(df,x=(feature$group),y=(name_leave),palette=palette_multicolor,color=feature$group)
		print("yes")
	}
	else if(feature$fill != 'no' && feature$err_bar=='yes'){
		p <- ggboxplot(df,x=(feature$group),y=(name_leave),palette=palette_multicolor,color=feature$group,add='jitter')
	}
	else{
		p <- ggboxplot(df,x=(feature$group),y=(name_leave),palette=palette_multicolor,add='jitter')
	}

	if(feature$markers == 'no'){
		return (p)
	}
	else{
		#p <- p + stat_compare_means()

		if(feature$asterisk == 'no'){
			p <- p+ stat_compare_means(comparisons=comparison,method=feature$smethod)
		}
		else{
			p <- p + stat_compare_means(comparisons=comparison,label='p.signif',method=feature$smethod)
		}

	}
}

barcharts_plot <- function(df,feature,comparison){
	library(ggsignif)
	nam_all <- names(df)
	name_leave = nam_all[which(nam_all!=feature$group)]
	name_leave = name_leave[which(name_leave!=feature$subgroup)]

	if(feature$markers == 'no'){
		

		
			p <- ggplot(data=df,aes(get(feature$group),get(name_leave[1]),fill=get(feature$subgroup))) + geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge())

		if(feature$err_bar != 'no'){
			#p <- p + stat_summary(fun.min=err_min,fun.max=err_max, geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))
			p <- p + stat_summary(fun.max=err_max, fun.min=mean,geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))
		}
		p <- p + labs(x=feature$group,y=name_leave,fill=feature$subgroup)
		return(p)
	}
	else{
		p <- ggplot(data=df,aes(get(feature$subgroup),get(name_leave[1])))
		if (feature$fill == 'no'){
			p <- p + stat_summary(geom="bar",fun=mean,color='black') 
		}
		else{
			p <- p + stat_summary(geom="bar",fun=mean,aes(fill=get(feature$fill)),color='black')
		}

		if(feature$err_bar != 'no'){
			#p <- p + stat_summary(fun.min=err_min,fun.max=err_max, geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))
			p <- p + stat_summary(fun.max=err_max, fun.min=mean,geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))
		}

		p <- p + facet_wrap(~get(feature$group),nrow=1,strip.position='bottom')
		if(feature$asterisk == 'no'){
			feature$asterisk=F
		}
		else{
			feature$asterisk=T
		}
		print('yes we can')
		print(str(comparison))
		print('am i')
		print(str(df))
		comparison <- lapply(comparison,as.numeric)
		p <- p + geom_signif(comparisons=comparison,step_increase=feature$interval,map_signif_level =feature$asterisk,size=feature$size,textsize=feature$textsize,test=feature$smethod) + theme(axis.text.x=element_blank(),panel.spacing.x=unit(0,'pt'),strip.background=element_blank(),axis.ticks.length =unit(0,'pt')) + labs(x=feature$group,y=name_leave,fill=feature$subgroup)
		return(p)
	}
}

boxcharts_plot <- function(df,feature,comparison){
	library(ggpubr)
	nam_all <- names(df)
	name_leave = nam_all[which(nam_all!=feature$group)]
	name_leave = name_leave[which(name_leave!=feature$subgroup)]
	if(feature$markers == 'no'){
		if(feature$fill == 'no' && feature$err_bar=='no'){
			

			p <- ggboxplot(df,x=(feature$group),y=name_leave,palette='jco')
		}
		else if(feature$fill != 'no' && feature$err_bar=='no'){
			p <- ggboxplot(df,x=(feature$group),y=name_leave,palette='jco',color=(feature$fill))
		}
		else if(feature$fill != 'no' && feature$err_bar=='yes'){
			p <- ggboxplot(df,x=(feature$group),y=name_leave,palette='jco',color=(feature$fill),add='jitter')
		}
		else{
			p <- ggboxplot(df,x=(feature$group),y=name_leave,palette='jco',add='jitter')
		}			
	}
	else{
		if(feature$fill == 'no' && feature$err_bar=='no'){
			p <- ggboxplot(df,x=(feature$subgroup),y=name_leave,palette='jco') + facet_wrap(~get(feature$group),nrow=1,strip.position='bottom')+ theme_grey() +theme(axis.text.x=element_blank(),panel.spacing.x=unit(0,'pt'),strip.background=element_blank(),axis.ticks.length =unit(0,'pt'))
			if(feature$asterisk == 'no'){
				p <- p+ stat_compare_means(comparisons=comparison,method=feature$smethod)
			}
			else{
				p <- p + stat_compare_means(comparisons=comparison,label='p.signif',method=feature$smethod)
			}
			
		}
		else if(feature$fill != 'no' && feature$err_bar=='no'){
			p <- ggboxplot(df,x=(feature$subgroup),y=name_leave,palette='jco',color=(feature$fill))+ facet_wrap(~get(feature$group),nrow=1,strip.position='bottom')+ theme_grey() +theme(axis.text.x=element_blank(),panel.spacing.x=unit(0,'pt'),strip.background=element_blank(),axis.ticks.length =unit(0,'pt'))
			if(feature$asterisk == 'no'){
				p <- p+ stat_compare_means(comparisons=comparison,method=feature$smethod)
			}
			else{
				p <- p + stat_compare_means(comparisons=comparison,label='p.signif',method=feature$smethod)
			}
		}
		else if(feature$fill != 'no' && feature$err_bar=='yes'){
			p <- ggboxplot(df,x=(feature$subgroup),y=name_leave,palette='jco',color=(feature$fill),add='jitter')+ facet_wrap(~get(feature$group),nrow=1,strip.position='bottom')+ theme_grey() +theme(axis.text.x=element_blank(),panel.spacing.x=unit(0,'pt'),strip.background=element_blank(),axis.ticks.length =unit(0,'pt'))
			#p <- ggboxplot(df,x=(feature$subgroup),y=name_leave,palette='jco',color=(feature$fill),add='jitter')
			if(feature$asterisk == 'no'){
				p <- p+ stat_compare_means(comparisons=comparison,method=feature$smethod)
			}
			else{
				p <- p + stat_compare_means(comparisons=comparison,label='p.signif',method=feature$smethod)
			}
		}
		else{
			p <- ggboxplot(df,x=(feature$subgroup),y=name_leave,palette='jco',add='jitter')+ facet_wrap(~get(feature$group),nrow=1,strip.position='bottom')+ theme_grey() +theme(axis.text.x=element_blank(),panel.spacing.x=unit(0,'pt'),strip.background=element_blank(),axis.ticks.length =unit(0,'pt'))
			if(feature$asterisk == 'no'){
				p <- p+ stat_compare_means(comparisons=comparison,method=feature$smethod)
			}
			else{
				p <- p + stat_compare_means(comparisons=comparison,label='p.signif',method=feature$smethod)
			}
		}			
	}	
}

construct_allcomp <- function(data,feature){
	group <- as.character(data[[feature$group]])
	uniq_group <- group[!duplicated(group)]
	if (length(uniq_group)==1){
		compaired=list()
		return (compaired)
	}
	if(length(feature$compare_first) != 0){
		compaired=list()
		for( i in 1:length(uniq_group)){
			if(i != 1){
				compaired <- c(compaired,list(c(uniq_group[1],uniq_group[i])))
				
			}
		}
		print(compaired)
		return(compaired)
	}
	compaired <- combn(uniq_group,2,simplify=F)
	return (compaired)
}

bar_mean <- function(data,feature,subtype){
	print('1')
	comparisons <- construct_allcomp(data,feature)
	print("2")
	p <- barmean_plot(data,feature,comparisons)
	print("3")
	p <- p + theme_light() + theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
	
	plot_continue(p,feature)
}

box_mean <- function(data,feature,subtype){
	comparisons <- construct_allcomp(data,feature)
	p <- boxmean_plot(data,feature,comparisons)
	plot_continue(p,feature)
}


construct_charts <- function(data,feature){
	#全部都做
	sub <- unique(data[[feature$subgroup]])
	print(str(sub))
	if(length(feature$compare_first) != 0){
		compaired=list()
		for( i in 1:length(sub)){
			if(i != 1){
				compaired <- c(compaired,list(c(sub[1],sub[i])))
				print(str(compaired))
			}
		}
		return(compaired)
	}


	compaired <- combn(sub,2,simplify=F)
	return(compaired)
}

bar_charts <- function(data,feature,subtype){
	comparisons <- construct_charts(data,feature)
	p <- barcharts_plot(data,feature,comparisons)
	plot_continue(p,feature)
}

box_charts <- function(data,feature,subtype){
	comparisons <- construct_charts(data,feature)
	p <- boxcharts_plot(data,feature,comparisons)
	plot_continue(p,feature)
}

pcoa_main <- function(data,feature,subtype){
	res <- pcoa_count(data,feature,subtype)
	download_save <- pca_save_file(res,feature,data)
	p <- plotpcoa(download_save,feature)
	#feature$maxdetail=0
	plot_continue(p,feature)
	#main_file = paste0(feature$out,"/ordination.png")
	#ggsave(main_file,plot=p)
	#save_plot <- paste0(feature$out,"/plotmain")
	#save(p,file=save_plot)
}

pca_save_file <- function(res,feature,data){
    res_file = paste0(feature$out,"/res.Data")
    save(res,file=res_file)

    sre_file = paste0(feature$out,'/sre.txt')
    write.table(res$values,sre_file,sep='\t')

    score_file = paste0(feature$out,'/score.txt')
    write.table(as.data.frame(res$vectors),score_file,sep='\t')

    xlab_new <- paste0("PC1(",round(res$values$Relative_eig[1]*100,2),"%)")
    ylab_new <- paste0("PC2(",round(res$values$Relative_eig[2]*100,2),"%)")
    lab_new <- c(xlab_new,ylab_new)

    download_save = res$vectors[,1:2]
    download_save = as.data.frame(download_save)
    colnames(download_save)=lab_new
    rownames(download_save)= rownames(data)
    download_save$text = rownames(data)


    if (length(feature$group) != 0){
        download_save$group = data[[feature$group]]
		download_save$group = as.factor(download_save$group)
    }

#temp
	if (length(feature$double) != 0){

		download_save <- rbind(download_save,download_save)
		print("double")
	}
#temp


    download_file = paste0(feature$out,"/out.txt")
    write.table(download_save,download_file,sep='\t')

	return(download_save)
}

plotpcoa <- function(df,feature){
	print(str(df))
	if (length(df$group)!=0 && length(feature$typed)!=0){
		p1 <- ggplot(data=df,aes(get(colnames(df)[1]),get(colnames(df)[2]),color=group)) + stat_ellipse(aes(colour=group),type=feature$typed,level = as.numeric(feature$level),show.legend=F) + geom_point(size=3) + scale_colour_brewer(palette="Set1") + labs(x=colnames(df)[1],y=colnames(df)[2])
	}
	if (length(df$group)!=0 && length(feature$typed)==0){
		p1 <- ggplot(data=df,aes(get(colnames(df)[1]),get(colnames(df)[2]),color=group)) + geom_point(size=3) + scale_colour_brewer(palette="Set1") + labs(x=colnames(df)[1],y=colnames(df)[2])
	}
	if(length(df$group)==0){
		p1 <- ggplot(data=df,aes(get(colnames(df)[1]),get(colnames(df)[2]))) + geom_point(size=3) + labs(x=colnames(df)[1],y=colnames(df)[2])
	}	
	return(p1)
}

exgroup <- function(data,feature){
	#print(feature)
	if (length(feature$group) != 0){
		data <- data[,which(colnames(data)!=feature$group)]
	}
	return (data)
}

standard <- function(data,feature){
	if (feature$stand == "yes"){
		data <- apply(data,2,scale,center=T,scale=T)
	}
	return (data)
}


pcoa_count <- function(data,feature,subtype){
	library(vegan)
	library(ape)
	#print(feature)
	data <- exgroup(data,feature)
	data <- standard(data,feature)

	print(feature$distance)
	dis <- vegdist(data,feature$distance,na.rm=T)
	res <- pcoa(dis)
	if(subtype=="11"){
		sre <- paste0(feature$out,"/sreplot.png")
		bio <- paste0(feature$out,"/biplot.png")

		x_num <- length(colnames(res$vectors))
	
		data_bar <- data.frame(x=colnames(res$vectors),y=res$values$Eigenvalues[1:x_num])
		data_bar$x <- factor(data_bar$x,levels=data_bar$x)
		sre_plot <- ggplot(data_bar,aes(x,y)) + geom_bar(stat="identity",fill="black") + scale_y_continuous(expand=c(0,0)) + theme_classic()
		ggsave(sre,plot=sre_plot)
		
		#print(str(res))
		#print(str(data))

		png(bio)
		bio_plot = biplot(res,data)
		print(bio_plot)
		dev.off()
		#ggsave(bio,plot=bio_plot)
	}

	return (res)
}

main()
