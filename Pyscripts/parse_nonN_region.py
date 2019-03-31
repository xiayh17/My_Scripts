from pyfaidx import Fasta
import sys
def parse_nonN(fasta):
	fasta = Fasta(fasta,as_raw=True)
	chroms = fasta.keys()
	for chrom in chroms:
		seq = fasta[chrom][:]
		start = 0
		end = 0
		mark = 0
		end_index = len(seq)-1
		for index,base in enumerate(seq):
			if base == 'N':
				if mark == 1:
					print("%s\t%d\t%d" % (chrom, start, index))
					mark = 0
			else:
				if mark == 0:
					mark = 1
					start = index
				if index == end_index:
					print("%s\t%d\t%d" % (chrom, start, index+1))

if __name__ == '__main__':
	parse_nonN(sys.argv[1])
