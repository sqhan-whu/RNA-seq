de_gene <- as.vector(read.delim('data\\252.txt',head=F, sep = '\t')[[1]])
gene_list <- read.delim('C:\\Users\\50246\\Desktop\\scratch\\GO\\gencode.v31.gene_length.txt', sep = '\t')
genes <- rep(0, nrow(gene_list))

names(genes) <- gene_list[,1]

genes[de_gene] <- 1

pwf <- nullp(DEgenes = genes, bias.data = gene_list$length, plot.fit = FALSE)

pvals <- goseq(pwf, 'hg38', 'ensGene')
pvals$FDR <- p.adjust(pvals$over_represented_pvalue, method = 'fdr')

#pvals <- pvals[c(1,2,8,4,5,6,7)] #第三列的p值就不要了

getGeneLists <- function(pwf, goterms, genome, ids){

gene2cat <- getgo(rownames(pwf), genome, ids)

cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),

unlist(gene2cat, use.names = FALSE))

out <- list()

for(term in goterms){

tmp <- pwf[cat2gene[[term]],]

tmp <- rownames(tmp[tmp$DEgenes > 0, ])

out[[term]] <- tmp

}

out

}

goterms <- pvals$category

goList <- getGeneLists(pwf, goterms, 'hg38', 'ensGene')

#默认将富集到GO的基因名称添加在最后一列，然后输出到本地

pvals$Ensembl_ID <- sapply(pvals$category, function(x) paste(goList[[x]], collapse = '; '))

write.table(pvals, 'goseq2.txt', sep = '\t', row.names = FALSE, quote = FALSE)
