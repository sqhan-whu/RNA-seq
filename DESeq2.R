#####
htseq-count -f bam -s no -r pos -a 10 -t exon -i gene_id -m intersection-nonempty Input-con-P1.Aligned.sortedByCoord.out.bam gencode.v31.annotation.gtf


##  combine count table
paste process1/process/Input-con-P1/Input-con-P1.htseq-count-id.tab 
process1/process/Input-con-P2/Input-con-P2.htseq-count-id.tab 
process1/process/Input-con-P3/Input-con-P3.htseq-count-id.tab 
process2/process/inputP1/inputP1.htseq-count-id.tab 
process2/process/inputP2/inputP2.htseq-count-id.tab 
process2/process/inputP3/inputP3.htseq-count-id.tab 
| awk '{printf $1 "\t";for(i=2;i<=NF;i+=2) printf $i"\t";printf "\n"}' > combine.txt


library('ggpubr')
library('DESeq2')
library('biomaRt')
library("curl")
library('stringr')

mycounts<-read.table("diff_gan.txt",head=T,sep='\t')
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
condition <- factor(c(rep("control",2),rep("treat",3)), levels = c("control","treat"))
colData <- data.frame(row.names=colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = lfcShrink(dds, coef=resultsNames(dds)[2], type = 'apeglm')
res = res[order(res$pvalue),]
head(res)
summary(res)

#write.csv(res,file="All_gan_results.csv")

mydata <- as.data.frame(res)
mydata$Condition=ifelse(mydata$log2FoldChange>=1 & mydata$padj<=0.05,"up",ifelse(mydata$log2FoldChange<=-1 & mydata$padj<=0.05,"down","normal"))
mydata$color=ifelse(mydata$log2FoldChange>=1 & mydata$padj<=0.05,"red",ifelse(mydata$log2FoldChange<=-1 & mydata$padj<=0.05,"green","gray"))
diff_gene_deseq2 <-subset(mydata,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

#write.csv(diff_gene_deseq2, file="diff_gene_deseq2.csv")

#search database
## searchDatasets(mart = ensembl, pattern = "hsa.*")

ensembl = useDataset("hsapiens_gene_ensembl",useMart("ensembl"))

rownames(diff_gene_deseq2) <- str_sub(rownames(diff_gene_deseq2), 1, 15)

my_ensembl_gene_id <-rownames(diff_gene_deseq2)

my_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = ensembl)

ensembl_gene_id<-rownames(diff_gene_deseq2)

diff_gene_deseq2<-cbind(ensembl_gene_id, diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
diff_gene_deseq2<- diff_gene_deseq2[,-2]


diff_name<-merge(diff_gene_deseq2,my_symbols,by="ensembl_gene_id")
diff_name= diff_name[order(diff_name$pvalue),]
write.csv(diff_name, file="gan_diff_gene_deseq2_name_annotation.csv")




#mydata$ensembl_gene_id<- str_sub(rownames(mydata), 1, 15)

mydata$ensembl_gene_id<- rownames(mydata)
name <-read.table("gene_id_name.txt",head=T,sep='\t')
total_name <-merge(mydata,name,by="ensembl_gene_id")
total_name = total_name[order(total_name$pvalue),]
write.csv(total_name, file="gan_total_gene.csv")




rownames(total_name) <-total_name$external_gene_name


pdf(file = "diff_gan.pdf",width=8,height=8)
p4 <- ggmaplot(total_name, main = expression("Lung positive" %->% "negtive"),
               fdr = 0.05, fc = 2, size = 0.6,
               palette = c("#B31B21", "#1465AC", "darkgray"),
               #genenames = as.vector(diff_express$name),
               legend = "top", top = 25, 
               select.top.method = "fc",
               font.label = c("bold", 11), label.rectangle = TRUE,
               font.legend = "bold",
               font.main = "bold",
               ggtheme = ggplot2::theme_minimal())
p4
dev.off()













###########################################################################################################

## 3 control vs 3 treat

library('DESeq2')
library('biomaRt')
library("curl")
library('stringr')

mycounts<-read.table("combine.txt",head=F,sep='\t')
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
mycounts<-mycounts[,-7]
condition <- factor(c(rep("control",3),rep("treat",3)), levels = c("control","treat"))
colData <- data.frame(row.names=colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)

res = lfcShrink(dds, coef=resultsNames(dds)[2], type = 'apeglm')
res = res[order(res$pvalue),]
head(res)
summary(res)

write.csv(res,file="All_results.csv")

diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
#write.csv(diff_gene_deseq2, file="diff_gene_deseq2.csv")

#search database
## searchDatasets(mart = ensembl, pattern = "hsa.*")

ensembl = useDataset("hsapiens_gene_ensembl",useMart("ensembl"))


rownames(diff_gene_deseq2) <- str_sub(rownames(diff_gene_deseq2), 1, 15)

my_ensembl_gene_id <-rownames(diff_gene_deseq2)

my_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = ensembl)

ensembl_gene_id<-rownames(diff_gene_deseq2)

diff_gene_deseq2<-cbind(ensembl_gene_id, diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
diff_name<-merge(diff_gene_deseq2,my_symbols,by="ensembl_gene_id")
diff_name= diff_name[order(diff_name$pvalue),]
write.csv(diff_name, file="diff_gene_deseq2_name_annotation.csv")


