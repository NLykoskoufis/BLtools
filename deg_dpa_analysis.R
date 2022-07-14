#!/usr/bin/env Rscript 

library(ggplot2)
library(ggpubr)
library(DESeq2)
library(ggthemes)
library(tidyr)
library(ggvenn)
source("/Users/groups/braung/data/nikos/marleneMetabolism/binaries/common/hodgeLab.R")

###########################################
###### VARIABLES TO SET FOR ANALYSIS ######
###########################################

INPUT_FILE <- ""
OUTPUT_PREFIX <- ""
ROW_SUM_THRESHOLD <- 10 ### variable for filtering out lowly expressed genes / peaks

### LOAD MATRIX ### 
df <- as.data.frame(data.table::fread(INPUT_FILE, sep="\t", header=TRUE,stringsAsFactors = FALSE,))
# rename id column
names(df)[4] <- "id"

#### RUN PCA ANALYSIS ####

chrToKeepNoSexChrom <- paste0("chr",c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9"))

TMP <- df[df$Chr %in% chrToKeepNoSexChrom ,]
rownames(TMP) <- TMP$id
TMP <- as.matrix(TMP[,-c(1:6)])

colData <- data.frame("samples"=colnames(TMP))
colData$condition <- NA
colData$condition[grepl("", colData$samples)] <- ""
colData$condition[grepl("", colData$samples)] <- ""


colData$condition <- as.factor(colData$condition)
dds1 <- DESeqDataSetFromMatrix(countData = round(TMP,0),
                               colData = colData,
                               design= ~ condition)

vsd1 <- vst(dds1, blind=FALSE)
PCA <- plotPCA( vsd1, intgroup = "condition2",returnData=TRUE)


ggPCA <- plotPCA(vsd1, intgroup="condition") + 
  theme_clean()+
  scale_colour_tableau()+
  ggtitle("ATAC peaks PCA analysis")+
  theme(panel.grid=element_blank(), legend.position="top",legend.background = element_rect(color="black"),legend.title=element_blank(),plot.title=element_text(hjust=0.5))

ggsave(ggPCA, filename=paste(OUTPUT_PREFIX,".PCAplot.pdf",sep=""),height=5, width=5.5, device="pdf")

#### PERFORM DIFFERENTIAL PEAK / GENE ACTIVITY ##### 

chrToKeep <- paste0("chr",c("1","10","11","12","13","14","15","16","17","18","19","2","3","4","5","6","7","8","9","MT","X","Y"))
TMP <- df[df$Chr %in% chrToKeep ,]
rownames(TMP) <- TMP$id

colData <- data.frame("samples" = colnames(TMP))
colData$condition <- NA
colData$condition[grepl("", colData$samples)] <- "" # You need to fill with the required conditions. You can duplicate this line as many times as needed depending on the amount of conditions you have. 
colData$condition[grepl("", colData$samples)] <- ""

colData$condition <- as.factor(colData$condition)
dds1 <- DESeqDataSetFromMatrix(countData = round(TMP,0),
                               colData   = colData,
                               desig     = ~ condition)


keep <- rowSums(counts(dds1)) > ROW_SUM_THRESHOLD
dds1 <- dds1[keep,]
dds <- DESeq(dds1)
res <- results(dds)
norm.counts <- as.data.frame(counts(dds, normalized=TRUE))
norm.counts$id <- rownames(norm.counts)
norm.counts <- merge(df[,c(1:6)],as.data.frame(norm.counts),by="id") 
norm.counts <- dplyr::select(norm.counts, c(names(norm.counts)[2:4],"id",names(norm.counts)[5:length(names(norm.counts))]))

data.table::fwrite(norm.counts, file=paste(OUTPUT_PREFIX,".raw.counts.bed",sep=""),col.names=T,row.names=F,sep="\t", quote=FALSE)

pdf(paste(OUTPUT_PREFIX, ".foldChangeHistogram.pdf",sep=""),height=4.5, width=4)
hist(res$log2FoldChange, main="log2FoldChange",breaks=50)
dev.off()

pdf(paste(OUTPUT_PREFIX,".foldChangeHistogram.pdf",sep=""),height=4.5, width=8)
par(mfrow=c(1,2))
hist(res$pvalue, main = "pvalue distribution")
hist(res$padj, main = "adjusted pvalue distribution")
dev.off()

res$id <- rownames(res)
res <- merge(df[,c(1:6)],as.data.frame(res),by = "id")
res <- dplyr::select(res, c(names(res)[2:4],"id",names(res)[5:length(names(res))]))
res[is.na(res$padj),]$padj <- 1 


res$sig <- res$padj <= 0.05
res$log2FC_1_5_abs <- res$sig & abs(res$log2FoldChange) >= log2(1.5)
res$direction <- NA
res$direction[which(res$sig & res$log2FoldChange >0)] <- "upregulated_in_XXXXX" # change the XXXXXX with the desired condition.
res$direction[which(res$sig & res$log2FoldChange <0)] <- "downregulated_in_XXXXXX" # change the XXXXX with the desired condition
res$directionStringent <- NA 
res$directionStringent[which(res$log2FC_1_5_abs & res$log2FoldChange >0)] <- "upregulated_in_XXXXXX" # change the XXXXX with the desired condition
res$directionStringent[which(res$log2FC_1_5_abs & res$log2FoldChange <0)] <- "downregulated_in_XXXXXX" # change the XXXXX with the desired condition



#### Plot MA plot ####
png(paste(OUTPUT_PREFIX,"maplot_dpas.png",sep=""),height=5,width=6,units = "in",res=500)
maplot(res, pthresh = 0.05) 
dev.off()

#### Volcano plot ####
gg <- ggplot(res, aes(x= log2FoldChange, y = -log10(padj), colour=sig))+
  geom_point()+
  ggrepel::geom_text_repel(data=res[which(abs(res$log2FoldChange)>1.5 & -log10(res$padj)>=30),], aes(x= log2FoldChange, y = -log10(padj), label=id))+
  geom_hline(yintercept= -log10(0.05), color="red2", lty=2)+
  geom_vline(xintercept=log2(1.5), color="green", lty=2)+
  geom_vline(xintercept=-log2(1.5), color="green", lty=2)+
  scale_colour_manual(values=c("TRUE"="black", "FALSE"="lightgrey"))+
  labs(x=expression(log[2] * " fold change"), y= expression(-log[10] * " adjusted p-value"), title="SVZ vs DG\ndifferential peak activity")+
  theme_clean(base_size=15)+
  theme(legend.position=0,
        plot.title = element_text(hjust=0.5, face="bold"))

ggsave(gg, filename=paste(OUTPUT_PREFIX,"volcanoPlot_peaks_deseq2.pdf",sep=""),height=6,width=6.5,device="pdf")

#### WRITE DESEQ2 RESULTS TO FILE ####
data.table::fwrite(res, file=paste(OUTPUT_PREFIX,".DESeq2_results.txt",sep=""),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,na="NA")

upregulated <- res[which(res$log2FC_1_5_abs & res$log2FoldChange > 0),c(1:4)]
downregulated <- res[which(res$log2FC_1_5_abs & res$log2FoldChange < 0),c(1:4)]


##### Write upregulated genes or peaks

write.table(upregulated, file=paste(OUTPUT_PREFIX,"upregulated_in_XXXX.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") # change the XXXXX with the desired condition
write.table(downregulated, file=paste(OUTPUT_PREFIX,"downregulated_in_XXXX.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


###### ONLY FOR ATAC or ChIP-seq
# Write the coordinates of the peaks into a file to load in GREAT later one
write.table(res[,c(1:4),], sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE, file = paste(OUTPUT_PREFIX,".GREATready.bed",sep="")) ####

### plotting up/downregulated genes/peaks ### 

up <- res[which(res$padj <= 0.05 & res$log2FoldChange > 0),]
down <- res[which(res$padj <= 0.05 & res$log2FoldChange < 0),] 
rownames(norm.counts) <- norm.counts$id




pdf(paste(OUTPUT_PREFIX,"upregulated.pdf",sep=""),height=12,width=12)
par(mfrow=c(3,3))
for (rrr in 1:nrow(up))
  
{
  peak <- up[rrr,]$id
  cat("Processed:",rrr,"-",length(up),"\b")
  x <- data.frame("exp"=as.data.frame(t(norm.counts[peak,-c(1:6)])))
  x$region <- "DG"
  x$region[grepl("SVZ",rownames(x))] <- "SVZ"
  names(x)[1] <- "exp"
  boxplot(x$exp ~ x$region, frame=FALSE, col=c("#999999", "#E69F00"),main=paste0(peak,"\npval=",format.pval(up[rrr,]$padj)," log2FC=",round(up[rrr,]$log2FoldChange,2)) , xlab="", ylab="Normalized counts",ylim=c(0,max(x$exp)))
}
dev.off()



pdf(paste(OUTPUT_PREFIX,"downregulated.pdf",sep=""),height=12,width=12)
par(mfrow=c(3,3))
for (rrr in 1:nrow(down))
  
{
  peak <- down[rrr,]$id
  cat("Processed:",rrr,"-",nrow(down),"\b")
  x <- data.frame("exp"=as.data.frame(t(norm.counts[peak,-c(1:6)])))
  x$region <- "DG"
  x$region[grepl("SVZ",rownames(x))] <- "SVZ"
  names(x)[1] <- "exp"
  boxplot(x$exp ~ x$region, frame=FALSE, col=c("#999999", "#E69F00"),main=paste0(peak,"\npval=",format.pval(down[rrr,]$padj)," log2FC=",round(down[rrr,]$log2FoldChange,2)) , xlab="", ylab="Normalized counts",ylim=c(0,max(x$exp)))
  
}
dev.off()



