## usage: python Counts_To_FPKM.py HTSeq.table.txt ~/project/project/00.DATABASE/hg38/Annotation/gencode.v31.gene_length.txt > out.fpkm.txt

from sys import argv
import pandas as pd
pd.set_option('display.max_rows',None)
pd.set_option('display.max_colwidth', 1000)
pd.set_option('display.max_columns', 200)
pd.set_option('display.width', 1000)

def get_count(counts_path):
	gene_list = dict()
	col_sum = dict()
	df = pd.read_table(argv[1],sep='\s+')
	header = list(df)
	
	for i in range(len(df.sum())):
		if i ==0:
			continue
		else:
			col_sum[i] = df.sum()[i]

	row = df.shape[0]
	for row in range(df.shape[0]):
		if df.iloc[row,0]:
			gene_list[df.iloc[row,0]] = []
			for col in range(df.shape[1]):
				if col == 0:
					continue
				else:
					gene_list[df.iloc[row,0]].append(df.iloc[row][col])

	return gene_list, col_sum, header

def get_gene_len(gene_len_path):
	with open(gene_len_path) as f:
		gen_len = dict()
		for line in f:
			line = line.strip().split('\t')
			gen_len[line[0]] = line[1]
	return gen_len

def get_fpkm(gene_list, col_sum, gen_len, header):
	headers = header
	print('\t'.join(headers))
	for k, v in gene_list.items():
		new_row = ''
		fpkm = 0
		if k in gen_len:
#			new_row += str(k) + '\t' + str(gen_len[k])
			new_row += str(k)
			p =1
			for num in v:
				fpkm = (int(num) * 1000000000/(int(col_sum[p]) * int(gen_len[k])))
				new_row += '\t' + str(fpkm)
				p +=1
		print(new_row)


gene_list, col_sum, header  = get_count(argv[1])
gen_len = get_gene_len(argv[2])
get_fpkm(gene_list, col_sum, gen_len, header)
