BiocManager::install("MASS")

library('edgeR')
library('limma')
library('locfit')
library('splines')
library('KernSmooth')
library('statmod')
library('MASS')
library('dplyr')


rawdata <- as.matrix(read.delim('diff_fei.txt', sep = '\t', row.names = 1))
group = factor(c(rep(1, 2), rep(2,4)))
y = DGEList(counts=rawdata, genes=rownames(rawdata),group=group)
#y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)
keep = rowSums(cpm(y)>1) >= 2
y = y[keep,]
delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)
dge <- estimateDisp(y, design, robust = TRUE)


####### Normalize data : method 1
dgelist_norm <- calcNormFactors(y, method = 'TMM')
plotMDS(dgelist_norm, col = rep(c('red', 'blue'),2,5), dim = c(1, 2))
####### Normalize data : method 2 (log2 counts-permillion)
rawdata_norm = log(as.data.frame(t(t(rawdata)/colSums(rawdata)*1000000)),2)
plotMDS(rawdata_norm, col = rep(c('red', 'blue'),2,5), dim = c(1, 2))


#negative binomial generalized log-linear model 拟合
fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
lrt <- glmLRT(fit)   #统计检验
write.csv(topTags(lrt, n = nrow(y$counts)), 'glmLRT.csv', quote = FALSE)  


#quasi-likelihood negative binomial generalized log-linear model 拟合
fit <- glmQLFit(dge, design, robust = TRUE)        #拟合模型
lrt <- glmQLFTest(fit)    #统计检验
 
topTags(lrt)
write.csv(topTags(lrt, n = nrow(y$counts)), 'glmQLFTest.csv', quote = FALSE)        #输出主要结果
 
dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
summary(dge_de)
 
plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)




#Exact Test
et = exactTest(y)
et_tags = topTags(et, n=200)$table %>% filter(PValue < 0.001)

#GLM
fit = glmFit(y, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, n = 200)$table %>% filter(PValue < 0.001)

#QL F-test
fit = glmQLFit(y, design)
contrast.design=c(1,-1)
qlf = glmQLFTest(fit, contrast=contrast.design)
qlf_tags = topTags(qlf, n=200)$table  %>% filter(PValue < 0.001)

top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))




targets <- as.matrix(read.delim('diff_fei.txt', sep = '\t', row.names = 1))
group = factor(c(rep(1, 2), rep(2,4)))
dgelist <- DGEList(counts = targets, group = group)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep,]

dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))
design <- model.matrix(~group)    #构建分组矩阵
dge <- estimateDisp(dgelist_norm, design, robust = TRUE) #估算离散值
 
plotBCV(dge) #作图查看

#negative binomial generalized log-linear model 拟合
fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
lrt <- glmLRT(fit)   #统计检验
 
topTags(lrt)
write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'glmLRT.csv', quote = FALSE)        #输出主要结果
 
dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)  #查看默认方法获得的差异基因
summary(dge_de)
 
plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)





###########################################################################
edgeR:

targets <- as.matrix(read.delim('diff_fei.txt', sep = '\t', row.names = 1))
group = factor(c(rep(1, 2), rep(2,4)))

deglist <- DGEList(counts = targets, group= group)

keep = rowSums(cpm(deglist)>1) >= 2

deglist = deglist[keep,]
deglist$samples$lib.size = colSums(deglist$counts)

delivery = group
data.frame(Sample=colnames(deglist),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(deglist)

deglist = estimateGLMCommonDisp(deglist, design, verbose=TRUE)
deglist = estimateGLMTrendedDisp(deglist,design)
deglist = estimateGLMTagwiseDisp(deglist,design)

#Exact Test
et = exactTest(deglist)
et_tags = topTags(et, nrow(dgelist$counts))$table %>% filter(PValue < 0.001)

#GLM
fit = glmFit(deglist, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, nrow(dgelist$counts))$table %>% filter(PValue < 0.001)

#QL F-test
fit = glmQLFit(deglist, design)
contrast.design=c(1,-1)
qlf = glmQLFTest(fit, contrast=contrast.design)
qlf_tags = topTags(qlf, nrow(dgelist$counts))$table  %>% filter(PValue < 0.001)

#Merge tests
top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))

rawdata_norm = as.data.frame(t(t(targets)/colSums(targets)*1000000))



dgelist_norm <- calcNormFactors(y, method = 'TMM')
plotMDS(dgelist_norm, col = rep(c('red', 'blue'),2,5), dim = c(1, 2))

design <- model.matrix(~group)
dge <- estimateDisp(dgelist_norm,design, robust=TRUE)
plotBCV(dge)

fit <-glmFit(dge,design,robust=TRUE)
lrt <-glmLRT(fit)
write.csv(topTags(lrt, n =nrow(dgelist$counts)),'glmLRT.csv',quote = FALSE)

dge_de <-decideTestsDGE(lrt,adjust.method='fdr', p.value=0.05)
summary(dge_de)

plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     #作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)

fit = glmQLFit(dge,design,robust=TRUE)
lrt = glmQLFTest(fit)
write.csv(topTags(lrt, n =nrow(dgelist$counts)),'glmQLFTest.csv',quote = FALSE)
dge_de <-decideTestsDGE(lrt,adjust.method='fdr', p.value=0.05)
summary(dge_de)







rawdata = read.csv("diff_fei.txt")
group = factor(c(rep(1, 2), rep(2,4)))  #1 = Term, 2 = Preterm

y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)
keep = rowSums(cpm(y)>1) >= 2
y = y[keep,]

y$samples$lib.size = colSums(y$counts)

delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)

y = estimateGLMCommonDisp(y, design, verbose=TRUE)
y = estimateGLMTrendedDisp(y,design)
y = estimateGLMTagwiseDisp(y,design)

#Exact Test
et = exactTest(y)
et_tags = topTags(et, n =nrow(y$counts))$table %>% filter(PValue < 0.01)

#GLM
fit = glmFit(y, design)
contrast.design= c(-1,1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, n =nrow(y$counts))$table %>% filter(PValue < 0.01)

#QL F-test
fit = glmQLFit(y, design)
contrast.design= c(-1,1)
qlf = glmQLFTest(fit, contrast=contrast.design)

qlf_tags = topTags(qlf, n =nrow(y$counts))$table  %>% filter(PValue < 0.01)
#Merge tests
top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))
length(top_genes)
#Normalize data
rawdata_norm = as.data.frame(t(t(rawdata[,-1])/colSums(rawdata[,-1])*1000000))



####################################################################################################################
library('edgeR')
library('limma')
library('locfit')
library('splines')
library('KernSmooth')
library('statmod')
library('MASS')
library('dplyr')

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Get raw data
rawdata = read.csv("raw_data/rnaseq_counts.csv")
group = factor(c(rep(1, 7), rep(2,8)))  #1 = Term, 2 = Preterm

y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)
keep = rowSums(cpm(y)>1) >= 2
y = y[keep,]

y$samples$lib.size = colSums(y$counts)

delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)

y = estimateGLMCommonDisp(y, design, verbose=TRUE)
y = estimateGLMTrendedDisp(y,design)
y = estimateGLMTagwiseDisp(y,design)

#Exact Test
et = exactTest(y)
et_tags = topTags(et, n=200)$table %>% filter(PValue < 0.001)

#GLM
fit = glmFit(y, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, n = 200)$table %>% filter(PValue < 0.001)

#QL F-test
fit = glmQLFit(y, design)
contrast.design=c(1,-1)
qlf = glmQLFTest(fit, contrast=contrast.design)
qlf_tags = topTags(qlf, n=200)$table  %>% filter(PValue < 0.001)

#Merge tests
top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))

#Normalize data
rawdata_norm = as.data.frame(t(t(rawdata[,-1])/colSums(rawdata[,-1])*1000000))

#Filter data and save only genes of interest
idx = match(top_genes, rawdata$external_gene_name)
data_filtered= cbind(rawdata$external_gene_name[idx], rawdata_norm[idx,])
names(data_filtered)[names(data_filtered)=='rawdata$external_gene_name[idx]'] = 'genes'
data_filtered = data_filtered[-match('Y_RNA', data_filtered$genes),] #artifact of all female cohort
write.csv(data_filtered, 'fig_3/a/data_combined_edgeR.csv', row.names = FALSE)
