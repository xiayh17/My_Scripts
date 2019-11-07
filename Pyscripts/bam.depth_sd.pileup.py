""" COVERAGE DEPTH by Module from cnvlib.fix
   - Depth by pysam
   - Depth Standard Deviation including all bases on each target
   - Fast processing by multiprocess (concurrent MODULE)
"""

import pysam
import sys
import os
import statistics
import atexit
import tempfile
import gzip
from concurrent import futures
import numpy as np
import pandas as pd


def region_depth_count(bamfile, chrom, start, end):
	'''
	output:  chr1	12034	12952	12035	27
		 chr1	12034	12952	12036	25
		 '''
	def filter_read(read):
		return not (read.is_duplicate
				or read.is_secondary
				or read.is_unmapped
				or read.is_qcfail
				)
				#or read.mapq < min_mapq)
	count = 0
	bases = 0
	depth_d = dict()
	start = int(start)
	end = int(end)
	chrom = str(chrom)
	cov_base,coverage = 0,0
	for pos in range(start+1, end+1):
		depth_d.setdefault(pos, 0)
		
	#for pileupcolumn in bamfile.pileup(reference=chrom, start=start, end=end, truncate=True, flag_filter=800, min_base_quality=10, min_mapping_quality=1, ignore_orphans=True):
	for pileupcolumn in bamfile.pileup(reference=chrom, start=start, end=end, truncate=True):
		pos = pileupcolumn.reference_pos + 1      #0 based
		if pos in depth_d:
			depth_d[pos] = pileupcolumn.nsegments
			cov_base+=1
	#depth = bases / (end - start) if end > start else 0
	#row = (chrom, start, end, depth)
	#return count, row
	coverage = float(cov_base/len(depth_d))
	row = (chrom, start, end ,depth_d,coverage)
	return row


def split_file(bed_fname, chunk_size=5000):
	"""Split a BED file into `chunk_size`-line parts for parallelization."""
	import os
	k, chunk = 0, 0
	dirname = os.getcwd()
	fd, name = tempfile.mkstemp(suffix=".bed", prefix="tmp.%s." % chunk, dir=dirname)
	outfile = os.fdopen(fd, "w")
	atexit.register(rm, name)
	opener = (gzip.open if bed_fname.endswith(".gz") else open)
	with opener(bed_fname) as infile:
		for line in infile:
			if line[0] == "#":
				continue
			k += 1
			outfile.write(line)
			if k % chunk_size == 0:
				outfile.close()
				yield name
				chunk += 1
				fd, name = tempfile.mkstemp(suffix=".bed",prefix="tmp.%s." % chunk, dir=dirname)
				outfile = os.fdopen(fd, "w")
	outfile.close()
	if k % chunk_size:
		outfile.close()
		yield name

def split_dataframe(df, split_n=2):
	import numpy as np
	split_df = np.array_split(df, split_n)
	return(split_df)


def bed_cov(df, bam_fname, procs=2):
	split_df = split_dataframe(df, procs)
	print("-------split OK---------------")
	chunks = []
	with futures.ProcessPoolExecutor(procs) as pool:
		args_iter = ((bam_fname, subr) for subr in split_df)   ##### this part is critical, otherwise futures will pop out an ERROR "Can't pickle pysam object..."
		#print(list(args_iter))                                ##### So here a bam file, not pysam object should be provided
		for chunk in pool.map(_rdc, args_iter):                ##### But a pandas.dataframe can be provided directly
			chunks.extend(chunk)
	return chunks

def _rdc(args):
	"""Wrapper for parallel."""
	return _rdc_chunk(*args)

def _rdc_chunk(bamfile, bed_df):
	chunk = []
	import pysam
	#bamfile = pysam.Samfile(bamfile, 'rb')
	bamfile = pysam.AlignmentFile(bamfile, 'rb')
	for index, row in bed_df.iterrows():
		chrom, start, end, gc = row['chrom'],row['start'],row['end'],row['gc']
		chrom, start, end ,depth_d,coverage = region_depth_count(bamfile,chrom, start, end)
		if len(depth_d) < 2:
			norm_sd = 10000.000

			try:
				meandepth = float(sum(list(depth_d.values()))/len(list(depth_d)))
			except:
				meandepth = 0
		else:
			meandepth = float(sum(list(depth_d.values()))/len(list(depth_d)))
			sd = statistics.stdev(list(depth_d.values()))
			if meandepth == 0:
				norm_sd = 1000.000
			else:
				norm_sd = float(sd/meandepth)
		chunk.append((chrom, start, end, gc, meandepth, norm_sd, coverage))
	return chunk
			
def rm(path):
	"""Safely remove a file."""
	try:
		os.unlink(path)
	except OSError:
		pass



if __name__== '__main__':
	bamfile = sys.argv[1]
	bedfile = sys.argv[2]
	outfile = sys.argv[3]
	procs = int(sys.argv[4])
	#bed_df = pd.read_csv(bedfile, sep="\t", names=['chrom','start', 'end', 'gc'])
	bed_df = pd.read_csv(bedfile, sep="\t", names=['chrom','start', 'end', 'gc'],low_memory=False)
	bed_df.chrom = bed_df.chrom.astype('object')
	chunks = bed_cov(bed_df, bamfile, procs)
	df = pd.DataFrame.from_records(chunks, columns=['chrom','start', 'end', 'gc', 'meandepth','sd','coverage'])
	df.to_csv(outfile, index=False, sep='\t')
