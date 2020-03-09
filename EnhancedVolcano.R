library('DESeq2')
library('biomaRt')
library("curl")
library('stringr')
library('EnhancedVolcano')

mycounts<-read.table("combine2.txt",head=F,sep='\t')
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

#write.csv(res,file="All_results.csv")


result <- as.data.frame(res)

result$Gene <- rownames(result)

write.table(result,file = "N061011_vs_N61311_difference.txt",sep = "\t",row.names = TRUE,quote =FALSE)

inputFileName <- "N061011_vs_N61311_difference.txt"

volcanodata <- read.table(inputFileName, header=TRUE, sep="\t",stringsAsFactors = FALSE)

p <- EnhancedVolcano(volcanodata,
                lab = rownames(volcanodata),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-4, 4),
                ylim = c(0, -log10(10e-17)),
                title = 'normal versus patients',
                #pCutoff = 10e-4,
                #FCcutoff = 1.333,
                cutoffLineCol = 'white',
                transcriptPointSize = 1.2,
                transcriptLabSize = 3,
                labhjust = 0,
                vline = 0,
                hline = NULL,
                shadeFill = 'white',
                gridlines.major=FALSE,
                gridlines.minor =FALSE,
                #caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                col = c("lightgrey", "#a6b1e1", "#5b8c85", "#ff5151"),
                colAlpha = 0.8,
                drawConnectors = FALSE,
                borderColour = "black")

ggsave(p,file="N061011_vs_N61311_EnhancedVolcano_plot.pdf",width = 9,height = 9)
