import sys

def get_gc(fasta,bed):
	from pyfaidx import Fasta
	chroms = Fasta(fasta,as_raw=True)

	with open(bed,'r') as fh:
		for line in fh:
			line = line.rstrip()
			if line:
				lines = line.split('\t')
				chrom,start,end = lines
				start = int(start)
				end = int(end)
				seq = chroms[chrom][start:end]
				seq = seq.upper()
				A_cnt,T_cnt,C_cnt,G_cnt,N_cnt = (seq.count('A'),seq.count('T'),seq.count('C'),seq.count('G'),seq.count('N'))
				#total = sum(A_cnt,T_cnt,C_cnt,G_cnt,N_cnt)
				#C_cnt,G_cnt = (seq.count('C'),seq.count('G'))
				C_G_cnt = C_cnt+G_cnt
				ATCG_cnt = A_cnt+T_cnt+C_cnt+G_cnt
				#print("%s\t%.6f" % (line, C_G_cnt/ATCG_cnt))
				print("%s\t%.4f" % (line, float(C_G_cnt/ATCG_cnt)*100))
if __name__ == '__main__':
	fasta = sys.argv[1]
	bed = sys.argv[2]
	get_gc(fasta,bed)
