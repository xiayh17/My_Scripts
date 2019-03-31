import pysam
import os,sys

def get_mate(primary_rec, bamfile):   ## primary
	pass

def get_dup(primary_rec, bamfile, out_handle):
	duplist = []
	prim_chrom = primary_rec.reference_name
	prim_start = primary_rec.reference_start
	prim_qseq = primary_rec.seq
	prim_name = primary_rec.qname
	
	prim_mate_rec = bamfile.mate(primary_rec)    ## read mate
	if prim_start > 0:
		dup_recs = bamfile.fetch(prim_chrom, start=prim_start-1, end=prim_start+1)
	else:
		dup_recs = bamfile.fetch(prim_chrom, start=prim_start, end=prim_start+1)
	dup_recs = list(dup_recs)                    ## possible dups
	if dup_recs:
		for rec in dup_recs:
			if rec.reference_start == prim_start and rec.is_duplicate and rec.seq == prim_qseq and not rec.mate_is_unmapped:
				dup_mate = bamfile.mate(rec)
				if dup_mate:
					if dup_mate.reference_start == prim_mate_rec.reference_start and dup_mate.seq == prim_mate_rec.seq:
						duplist.append(rec)
	#print(prim_name)
	if duplist:
		name_list = [rec.qname for rec in duplist]
		#print(prim_name, '--'.join(name_list))
		out_handle.write(prim_name + ' ' + '--'.join(name_list) + '\n')
		out_handle.flush()


if __name__== '__main__':
	bam = sys.argv[1]
	out = sys.argv[2]
	bamfile = pysam.Samfile(bam, 'rb')
	out_handle = open(out, 'w')

	for read in bamfile.fetch():
		#not read.is_secondary not read.is_unmapped
		if read.is_duplicate or read.is_supplementary or read.is_unmapped or read.is_secondary or read.is_read2 or read.mate_is_unmapped:
			pass
		else:
			#print('hahha----',read.qname)
			get_dup(read, bamfile, out_handle)
	out_handle.close()
