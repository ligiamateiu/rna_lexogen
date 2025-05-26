rm(list=ls(all=TRUE)) 
library(yaml)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(Hmisc)
library(NOISeq)
library(magrittr)

# library(pathview)
# ====================== define the function of DEA ======================

generate_dds<-function(inpath, outpath){
filelist<-paste(inpath,"/",samples.all,"_count.tsv",sep="")
gene.ids <- read.table(filelist[1], header=FALSE, sep="\t")[,1] 
count.table  <- do.call(cbind,lapply(filelist,function(fn)suppressWarnings(as.numeric(read.table(fn,header=FALSE, sep="\t")[,2]))))
row.names(count.table)<-gene.ids
colnames(count.table)<-samples.all


write.table(merge(count.table,gtf_symbols, by=0), file = paste(outpath,'/table_input_gene_counts.tsv',sep=""), quote=F,col.names=NA,row.names=T,sep="\t")


if(exclude.genes!="NO"){
	count.table<-count.table[rownames(count.table) %nin% exclude.genes,]
	print (paste("genes removed from the analysis:",exclude.genes, sep="\t"))
}
 
  # Filter protein-coding genes
  initnrgenes <- nrow(count.table)
  print(paste("Initial number of genes", initnrgenes, sep = " "))
  if (protein_coding_only == "YES") {
    ProteinCodingGenes <- gtf[gtf$gene_biotype == "protein_coding", ]
    count.table <- count.table[rownames(count.table) %in% ProteinCodingGenes$gene_id, ]
    print(paste("Kept only the protein coding genes:", nrow(count.table), sep = " "))
  }

  # Remove rRNA genes
  ribosomal_genes <- gtf[grepl("rRNA", gtf$gene_name, ignore.case = TRUE), ]
  count.table <- count.table[!rownames(count.table) %in% ribosomal_genes$gene_id, ]
  print(paste0("rRNA removed:", length(ribosomal_genes), "  Genes remaining:", nrow(count.table)))

  # Remove mtRNA genes
  mtRNA<-  gtf[seqnames(gtf) %in% c("MT", "M", "mt")]
  count.table <- count.table[!rownames(count.table) %in% mtRNA$gene_id, ]
  print(paste0("mtRNA removed:", length(mtRNA), "  Genes remaining:", nrow(count.table)))
  
  
    # Remove Y chromosome genes
  yRNA <- gtf[seqnames(gtf) %in% "Y"]
  count.table <- count.table[!rownames(count.table) %in% yRNA$gene_id, ]
  print(paste0("chrY RNA removed:", length(yRNA), "  Genes remaining:", nrow(count.table)))



    ## create the DESeqDataSet
    if(length(covariate1)>0 && length(covariate2)>0){
		 colData = data.frame(samples.all, group.all,covariate1.all,covariate2.all)
		dds.obj <- DESeqDataSetFromMatrix(countData = count.table, colData = colData, design = ~covariate2.all+covariate1.all+group.all)
	} else if (length(covariate1)>0 && length(covariate2)==0){
		colData = data.frame(samples.all, group.all,covariate1.all)
		dds.obj <- DESeqDataSetFromMatrix(countData = count.table, colData = colData, design = ~covariate1.all+group.all)
	} else{
		colData = data.frame(samples.all, group.all)
		dds.obj <- DESeqDataSetFromMatrix(countData = count.table, colData = colData, design = ~group.all)
	}
	
    #visualize sizefactors
    dds.obj <- estimateSizeFactors(dds.obj)
	size_factors <- sizeFactors(dds.obj)
    print(size_factors)
	size_factors_df <- data.frame(Sample = names(size_factors), SizeFactor = size_factors)
	p <- ggplot(size_factors_df, aes(x = Sample, y = SizeFactor)) + geom_bar(stat = "identity",fill = "blue") + labs(title = "Size Factors Across Samples", x = "Sample", y = "Size Factor") + geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  geom_hline(yintercept = 1.5, linetype = "dashed", color = "green") +   labs(title = "Size Factors Across Samples", x = "Sample", y = "Size Factor") + theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 6),  panel.background = element_rect(fill = "white"))
	print(p)
	ggsave(paste(outpath, "/plots/size_factors_plot.png",sep=''), plot = p, width = 8, height = 6, dpi = 300)
	

    
	keep <- rowSums( counts(dds.obj, normalized=TRUE) >= 10 ) >= 3
	
    normalized_counts <- counts(dds.obj, normalized=TRUE)
    write.table(normalized_counts,paste(outpath, '/library_normalized_gene_counts.tsv', sep = ''),quote=F,col.names=NA,row.names=T,sep="\t")


    dds.obj <- dds.obj[keep,]

    ## perform DEA
    dds.obj <- DESeq(dds.obj)
    save(dds.obj,file=paste(output.path,'/deasession.RData',sep=""))
	dim(dds.obj)
  return(dds.obj)
  
}

generate_allsamples_plots<-function(dds.obj){
    vst.obj <- vst(dds.obj, blind=F)
  
  pca_data <- assay(vst.obj)
  pca <- prcomp(t(pca_data))  # Transpose so samples are rows
  percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- colnames(pca_data)  # Add sample names
pca_df$Condition <- colData(dds.obj)$group.all

 pca_plot<- ggplot(pca_df, aes(x = PC1, y = PC2, color = group.all)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot") +
   theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    panel.grid.major = element_line(color = "gray90"),  # Light grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
ggsave(paste(output.path,'/plots/PCA_plot.png',sep=''), plot = pca_plot, width = 7, height = 7, dpi = 300)

#plotting samples 
  png(paste(output.path,'/plots/','SamplesDist_Euclidean.png',sep=''),width = 3600, height = 2400,res = 300)
  sampleDists <- dist(t(assay(vst.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vst.obj$samples.all, vst.obj$group.all, sep="_")
  colnames(sampleDistMatrix) <- paste(vst.obj$samples.all, vst.obj$group.all, sep="_")
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
  dev.off()
    return(vst.obj)

}


DEA <- function(control, treat,dds.obj,vst.obj) 	{
    # export the results

    res.dea <- results( dds.obj, contrast = c("group.all",control,treat))
    res.dea <- res.dea[complete.cases(res.dea), ]  # remove any rows with NA
# 
# 	ihw_res <- ihw(pvalue ~ baseMean,  data = as.data.frame(res.dea), alpha = fdr.thr)
#     rejections(ihw_res)
#     c(nbins(ihw_res), nfolds(ihw_res))
#     res.dea$ihw=adj_pvalues(ihw_res)
    
    dea <- as.data.frame(res.dea)
    dea <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    table_logs<-data.frame(assay(vst.obj),check.names=FALSE)
    table_logs <- merge(table_logs, gtf_symbols, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
	write.table(table_logs,paste(output.path, '/vst_normalized_gene_counts.tsv',sep=''),quote=F,col.names=TRUE,row.names=FALSE,sep="\t")

	
	deg <- dea[dea$padj < fdr.thr, ]
    if (nrow(deg) > 1) {
      deg <- deg[order(deg$padj, -abs(deg$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    }

    # save the DEAs and DEGs to files
	deg <- merge(deg, gtf_symbols, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
    dea <- merge(dea, gtf_symbols, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
    
    
    write.table(dea, paste(output.path, '/LogRatios_', control, '_', treat, '.tsv', sep = ''), quote=F,col.names=TRUE,row.names=FALSE, sep = '\t')
    write.table(deg, paste(output.path, '/SignifDiffExpr_', control, '_', treat, '.tsv', sep = ''), quote=F,col.names=TRUE,row.names=FALSE, sep = '\t')
    top20expr<-dea[order(dea$baseMean, decreasing = TRUE), ] 
    top20expr<-top20expr[1:20,]
	write.table(top20expr,paste(output.path, '/top20most_expressed_genes',control, '_', treat,'.tsv',sep=''),quote=F,col.names=TRUE,row.names=FALSE,sep="\t")

# plotting 
# plottingDEG 
# volcano plt 	
if(nrow(deg)>1){
	deg[["id"]] <- deg$gene_name
	show <- as.data.frame(deg[1:10, c("log2FoldChange", "padj", "id")])	
	deg_plot<-deg[!is.na(deg$id),]
	png(paste(output.path,"/plots/DEG_volcanoplot_enhanced_",control, "_", treat,".png",sep=""), width = 12, height = 8, units = 'in',res=500)
	p<-EnhancedVolcano(deg_plot,lab = deg_plot$gene_name,x = 'log2FoldChange',y = 'padj',pCutoff = fdr.thr,FCcutoff = fc.thr,cutoffLineType = 'twodash', cutoffLineWidth = 0.8,pointSize = 4.0,labSize = 6.0,legendLabels=c('Not sig.','log2FC','padj', 'padj & logFC'),legendPosition = 'right',legendLabSize = 16, legendIconSize = 5.0,title= paste("Volcano plot of signif DE genes at  BH adjusted p-value < ",fdr.thr,"log2FC > |",fc.thr,"|", sep=""),boxedLabels = TRUE,labCol = 'black',colAlpha = 4/5,drawConnectors = TRUE,widthConnectors = 1.0,colConnectors = 'black')
	plot(p)
	dev.off()
}


vsd_de<- data.frame(assay(vst.obj[rownames(vst.obj) %in% deg$Row.names,] ),check.names=FALSE)
#plot heatmap for the DEG if they are >=3
if(nrow(vsd_de)>2){
png(paste(output.path,"/plots/DEG_heatmap_",control, "_", treat,".png",sep=""),width = 2400, height = 1600,res = 300)
   vsd_de_plot<-vsd_de
   colsymbol<-dea[row.names(vsd_de_plot),which( colnames(dea)=="gene_name" )]
   rownames(vsd_de_plot)<-paste0(colsymbol,"_",rownames(vsd_de_plot))    
    df <- as.data.frame(colData(dds.obj)[,c("group.all")])
    colnames(df)<-"group"
    row.names(df)<-row.names(colData(dds.obj))
    pheatmap(vsd_de_plot, cluster_rows=TRUE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df)
dev.off()  
}

if (length(plot.genes.extra)>1) {
    dirplotextra <- paste(output.path, "/plots/plots_extra", sep = "")
    dir.create(dirplotextra, showWarnings = FALSE, recursive = TRUE)  # Create the directory if it doesn't exist

    normalized_counts <- counts(dds.obj, normalized=TRUE)
    ncounts_extra<-data.frame(normalized_counts[rownames(normalized_counts) %in% plot.genes.extra, ], check.names = FALSE)

    for (i in 1:ncol(ncounts_extra)) {
        extra_gene <- data.frame(group = colData(dds.obj)$group.all, gene = ncounts_extra[,i])
        ensid <- plot.genes.extra[i]
        colnames(extra_gene) <- c("group", "ncounts")
        
        fc <- round(dea[grep(ensid,dea[,1]), "log2FoldChange"], digits = 2)
        padj <- round(dea[grep(ensid,dea[,1]), "padj"], digits = 5)
        genesymbol <- dea[grep(ensid,dea[,1]), "gene_name"]
        
        if (length(genesymbol) == 0) {
            genesymbol <- "NA"
        }

        png(paste(dirplotextra, "/", genesymbol, "_", ensid, ".png", sep = ""), width = 12, height = 8, units = 'in', res = 500)
        p <- ggplot(extra_gene, aes(x = group, y = ncounts, fill = group)) +
            ggtitle(paste(genesymbol, ensid)) +
            geom_violin(trim = FALSE) +
            geom_boxplot(width = 0.1) +
            ylab("ncounts (library normalized counts)") +
            theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
                  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray"),
                  panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white"),
                  axis.text = element_text(size = 22),
                  axis.title = element_text(size = 22),
                  plot.title = element_text(size = 14),
                  legend.position = "none")
        print(p)
        dev.off()
    }
}

if(nrow(vsd_de)>0){
	dirplotdeg=paste(output.path,"/plots/plots_deg",sep="")
	dir.create(dirplotdeg)

##plot top 50 most differentially expressed genes

if(plotalldeg=="TOP"){
	topDEG<-deg$Row.names
	normalized_counts <- counts(dds.obj, normalized=TRUE)
	ncount_de_top<-data.frame(normalized_counts[rownames(normalized_counts) %in% topDEG, ], check.names = FALSE)

}
else{
	ncount_de_top <- counts(dds.obj, normalized=TRUE)
}


for(i in 1:nrow(ncount_de_top)){
	ncount_de_gene<-data.frame(group=colData(dds.obj)$group.all,gene=log(t(ncount_de_top[i,])))
	ensid<-rownames(ncount_de_top[i,])
	colnames(ncount_de_gene)<-c("group","ncount")	
	genesymbol<-deg[grep(ensid,deg[,1]),"gene_name" ]

    if(length(genesymbol)==0){genesymbol="NA"}
	png(paste(dirplotdeg,"/",genesymbol,"_",ensid,".png",sep=""), width = 12, height = 8, units = 'in',res=500)
 	p<-ggplot(ncount_de_gene,aes(x=group, y=ncount,fill=group)) + ggtitle(paste(genesymbol,ensid,sep="\t")) +geom_violin(trim=FALSE)+geom_boxplot(width=0.1)+ylab("log (normalized counts)")+theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"))+theme(axis.text=element_text(size=22),axis.title=element_text(size=22),plot.title = element_text(size = 14))+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

	print(p)
	dev.off()
    }
}
#covariates effect NOT YET
#resCov <- degCovariates(log2(counts(dds)+0.5),colData(dds))

  output<-dea
  return(output)

  }


# load the config file
yaml.file <- yaml.load_file('configs/config_main.yaml')

# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name of this analysis
controls <- yaml.file$CONTROL  # all groups used as control

covariate1<-yaml.file$COVARIATE1
covariate2<-yaml.file$COVARIATE2
print(controls)
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
print(treats)
fdr.thr <- yaml.file$FDR  # threshold of FDR/adjusted P-value for significantlly differentially expressed genes
meta.file <- yaml.file$METAFILE
organism<-yaml.file$ORGANISM
protein_coding_only<-yaml.file$PROTEIN_CODING_ONLY
fc.thr<-yaml.file$FCTHR
plotalldeg<-yaml.file$PLOTALLDEG
annotation<-yaml.file$ANNOTATION
plot.genes.extra<-yaml.file$PLOT_GENES
exclude.genes<-yaml.file$EXCLUDE_GENES_FROM_DEA
output.path <- file.path(yaml.file$FINALOUTPUT, "/dea")
input.path<-file.path(yaml.file$FINALOUTPUT, "/countFile")
# extract the metadata
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
meta.data <- meta.data[ which(meta.data$deause=='yes'), ]

group.all <- meta.data$group
samples.all<- meta.data$sample
samples.all

# Load GTF file
  gtf <- rtracklayer::import(annotation)
  
  # Extract gene symbols (assuming 'gene_name' is the attribute in the GTF)
 gtf_symbols <- unique(as.data.frame(gtf)[, c("gene_id", "gene_name")])
 row.names(gtf_symbols) <- gtf_symbols$gene_id


print(covariate1)
if(length(covariate1)>0){
covariate1.all<- meta.data[,covariate1]
}
if(length(covariate2)>0){
covariate2.all<- meta.data[,covariate2]
}

num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control




if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control
    
dds<-generate_dds(input.path,output.path)
vst<-generate_allsamples_plots(dds)

for (ith.comparison in c(1:num.comparison)) {
  control.set <- controls[ith.comparison]
  treat.set <- treats[ith.comparison]

### Do DEA
# the main function
data_fc<-DEA(control.set, treat.set,dds,vst)

    
}

sessionInfo()
