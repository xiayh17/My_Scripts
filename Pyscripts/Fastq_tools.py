import os,sys,re
import gzip
import io
from Bio.SeqIO.QualityIO import FastqGeneralIterator

### 获取反向互补序列
def trans_sequence(seq):
	'''get reversed Complementary sequence
	'''
	maps = str.maketrans('ATCG','TAGC')
	seqt = seq.translate(maps)
	seql = list(seqt)
	seql.reverse()
	seqs = ''.join(seql)
	return seqs

### 读取fastq 或 fastq.gz 文件

def read_fastq_Bio(fq):
	'''read fq or fq.gz file by FastqGeneralIterator
	   return generator
	'''
	if re.search(r'\.gz$',fqname):
		ft = gzip.open(fqname,'rb')
		fh = io.TextIOWrapper(ft)
	else:
		fh = open(fqname, 'r')
	
	fh_r = FastqGeneralIterator(fh)
	return fh_r
	fh.close()

def read_fastq(fq):
	'''read fq or fq.gz file by normal reader
	   return generator by yield
	'''
	if re.search(r'\.gz$',fqname):
		ft = gzip.open(fqname,'rb')
		fh = io.TextIOWrapper(ft)
	else:
		fh = open(fqname, 'r')
	
	while True:
		l1 = fh.readline()
		if not l1:
			break
		l2 = fh.readline()
		l3 = fh.readline()
		l4 = fh.readline()
		yield [l1,l2,l3,l4]

def filter_fastq(fqin,fqout):
	import gzip
	from io import StringIO
	from Bio import SeqIO,bgzf

	fo = gzip.open(fqout, "wt")
	with gzip.open(fqin, "rt") as handle:
		for record in SeqIO.parse(handle, "fastq"):
			# if record.id / record.seq / 
			fo.write(record.format("fastq"))

