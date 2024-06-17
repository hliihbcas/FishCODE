print(Sys.time())
library(methylKit)
library(genomation)
args=commandArgs(T)
inputfilenames = args[1]
samplenames = args[2]
mincov=args[3]

win.size=args[4]
step.size=args[5]
qvalue=args[6]
difference=args[7]
type=args[8]

treatment=args[10]

outfilename = args[9]
base_mincov = args[11]

bioprojects = args[12]

bathfla=args[13]
pca_per=args[14]
pvalue_pca=args[15]
batchids = args[16]
print(batchids)
batchids_c = strsplit(batchids,',')[[1]]
print(batchids_c)
treatment = as.numeric(strsplit(treatment,',')[[1]])
print(treatment)
print(str(treatment))

file.list = as.list(strsplit(inputfilenames,',')[[1]])

i = 1
for(x in file.list){
	
	temp = paste0('./data/example_meth/',x,'.methratio.txt')
	if(file.exists(temp)){
		file.list[[i]] = temp
		
	}else{
		temp = paste0('./data/example_meth/',x,'/common.data')
		file.list[[i]] = temp
	}		
	i=i+1
}

bed_file_name = paste0(outfilename,'.bed')

print(file.list)
samplenames.list = as.list(strsplit(samplenames,',')[[1]])

print(samplenames.list)
myobj <- methRead(
	file.list,		
	sample.id=samplenames.list,		
	assembly="own_design",
	treatment=treatment,
	context="CpG",
	mincov = as.numeric(mincov)
)
#myobj = filterByCoverage(myobj,lo.count=4)
print(Sys.time())
#print(attr(myobj,'class'))
myobj <- filterByCoverage(myobj,o.count=NULL,lo.perc=NULL,hi.count=NULL,,hi.perc=99.9)
#print(attr(myobj,'class'))
myobj <- normalizeCoverage(myobj, method = "median")
#print(attr(myobj,'class'))
#return(1)
print(as.numeric(mincov))	
print(as.numeric(step.size))	
print(as.numeric(win.size))
print(as.numeric(base_mincov))	
#	tileMethylCounts

sampleAnnotation = data.frame(batch_id=batchids_c)



regions = tileMethylCounts(myobj,win.size=as.numeric(win.size),step.size=as.numeric(step.size),cov.bases = as.numeric(base_mincov),mc.cores=16)

#print(attr(regions,'class'))
	

#print(head(regions[[1]]))
#return(1)	
meth = unite(regions,destrand = FALSE,mc.cores=16)
print(attr(meth,'class'))
if(bathfla=='batch' & length(unique(batchids_c))!=1){
	print(str(meth));print(str(sampleAnnotation));print(attr(meth,'class'))
	as = assocComp(mBase = meth,sampleAnnotation)
	pca_meth = data.frame(x=as$vars,y=as$association[1,])
	rownames(pca_meth)=seq(nrow(pca_meth))


	pcs= as.numeric(rownames(subset(pca_meth,x>as.numeric(pca_per) & y<as.numeric(pvalue_pca))))
	if(length(pcs)!=0){
		print("The following principal components met your deletion criteria and were deleted")
		print(paste('pc',as.character(pcs),sep=""))
		meth=removeComp(meth,comp=pcs)
		
	}else{
		print("*None of the principal components have a relationship with the batch effect that matches your deletion*")
	}
}

print(head(meth))	
if(bathfla=="covar" & length(unique(batchids_c))!=1){
	print("We treat data from different sources as covariates according to your choice. See calculateDiffMeth::covariates parameter for details")
	myDiff = calculateDiffMeth(meth,mc.cores=16,covariates=sampleAnnotation)
}else{

	myDiff = calculateDiffMeth(meth,mc.cores=16)
}

if (length(unique(batchids_c))==1){
	print("The samples you selected come from the same data set, so we did not detect and remove the batch effect")
	#myDiff = calculateDiffMeth(meth,mc.cores=4)
}

if (bathfla=="no"){
	print("We will not detect and remove batch effects based on your choice of no option")
	#myDiff = calculateDiffMeth(meth,mc.cores=4)
}


#myDiff = calculateDiffMeth(meth,mc.cores=8)
all = getMethylDiff(myDiff,difference=as.numeric(difference),qvalue=as.numeric(qvalue),type=type)
#gene.obj=readTranscriptFeatures("/home/hli/meth/test.bed",up.flank=2000,down.flank=0)
#diffAnn=annotateWithGeneParts(as(all,"GRanges"),gene.obj)
bedgraph(all,col.name = "meth.diff", file.name = bed_file_name)	
names(all)[3]='stop'


write.table(all,outfilename,sep="\t",quote=F,row.names=F)

#print(head(all))
#print(Sys.time())	
