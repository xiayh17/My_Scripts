from pyfaidx import Fasta
import sys

def parse_gaps(depth, Low_threshold, out_f):
	fi = open(out_f, 'w')
	with open(depth, 'r') as fh:
		chrom_s = 'NewChr1'
		mark=0
		pos_e = 0
		for line in fh:
			chrom,pos,dep = line.rstrip().split('\t')
			dep = int(dep)
			if not pos_e:
				pos_e = pos
			if chrom == chrom_s:
				if int(pos) - int(pos_e) > 1:
					if mark > 10:
						fi.write("%s\t%s\t%s\n" % (chrom_s, pos_s, pos_e))
						fi.flush()
					mark=0
				if dep < Low_threshold:
					if mark==0:
						chrom_s,pos_s=chrom,pos
						chrom_e,pos_e=chrom,pos
						mark=1
					else:
						chrom_e,pos_e = chrom,pos
						mark+=1
				else:
					if mark > 10:
						fi.write("%s\t%s\t%s\n" % (chrom_s, pos_s, pos_e))
						fi.flush()
					mark=0
			else:
				if mark > 10:
					fi.write("%s\t%s\t%s\n" % (chrom_s, pos_s, pos_e))
					fi.flush()
					mark=0
				if dep < Low_threshold:
					mark=1
					chrom_s,pos_s=chrom,pos
					chrom_e,pos_e=chrom,pos
				else:
					mark=0
				chrom_s = chrom
		if mark > 10:
			fi.write("%s\t%s\t%s\n" % (chrom_s, pos_s, pos_e))
			fi.flush()

if __name__ == '__main__':
	parse_gaps(sys.argv[1], int(sys.argv[2]), sys.argv[3])
