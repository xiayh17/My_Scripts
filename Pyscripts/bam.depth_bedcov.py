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
from skgenome.tabio import bedio


def main(bed_fname, bam_fname, by_count, min_mapq, processes):
  if by_count:
    #TODO
  elif by_bedcov:
    coverage_bedcov(bed_fname, bam_fname, by_count, min_mapq, processes)

      

 
def coverage_count(bed_fname, bam_fname, min_mapq, procs=1)
  regions = bedio.read(bed_fname)
  if procs == 1:
    bamfile = pysam.Samfile(bam_fname, 'rb')
    for chrom, subregions in by_chromosome(regions):
      logging.info("Processing chromosome %s of %s"
                   chrom, os.path.basename(bam_fname))
      for count, row in _rdc_chunk(bamfile, subregions, min_mapq):
        yield row
  else:
    with futures.ProcessPoolExecutor(procs) as pool:
      args_iter = ((bam_fname, subr, min_mapq) for _c, subr in regions.by_chromosome())
      for chunk in pool.map(_rdc, args_iter):
        for count, row in chunk:
          yield row
          
def _rdc(args):
    """Wrapper for parallel."""
    return list(_rdc_chunk(*args))
  
def _rdc_chunk(bamfile, regions, min_mapq):
    if isinstance(bamfile, basestring):
        bamfile = pysam.Samfile(bamfile, 'rb')
    for chrom, start, end, gene in regions.coords(["gene"]):
        yield region_depth_count(bamfile, chrom, start, end, gene, min_mapq)
 
def region_depth_fetch(bamfile, chrom, start, end, gene, min_mapq):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    def filter_read(read):
        """True if the given read should be counted towards coverage."""
        return not (read.is_duplicate
                    or read.is_secondary
                    or read.is_unmapped
                    or read.is_qcfail
                    or read.mapq < min_mapq)

    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            rlen = read.query_alignment_length
            if read.pos < start:
                rlen -= start - read.pos
            if read.pos + read.query_alignment_length > end:
                rlen -= read.pos + read.query_alignment_length - end
            bases += rlen
    depth = bases / (end - start) if end > start else 0
    row = (chrom, start, end, gene,
           math.log(depth, 2) if depth else NULL_LOG2_COVERAGE,
           depth)
    #return count, row
    return row

  def region_depth_pileup(bamfile, chrom, start, end, gene, min_mapq):
        """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    depth_d = dict()
    start = int(start)
	  end = int(end)
	  chrom = str(chrom)
	  cov_base,coverage = 0,0
	  for pos in range(start+1, end+1):
		  depth_d.setdefault(pos, 0)
		
	#for pileupcolumn in bamfile.pileup(reference=chrom, start=start, end=end, truncate=True, flag_filter=800, min_base_quality=10, min_mapping_quality=1,  
  for pileupcolumn in bamfile.pileup(reference=chrom, start=start, end=end,         ##defalt filter [UNMAP,SECONDARY,QCFAIL,DUP]
                                     truncate=True, min_mapping_quality=min_mapq,):
    pos = pileupcolumn.reference_pos + 1      #0 based
	  if pos in depth_d:
		  depth_d[pos] = pileupcolumn.nsegments
			#cov_base+=1
      
	#coverage = float(cov_base/len(depth_d))
	#row = (chrom, start, end ,depth_d,coverage)
  depth = float(sum(list(depth_d.values()))/len(list(depth_d)))
  row = (chrom, start, end, gene, math.log(depth, 2) if depth else NULL_LOG2_COVERAGE, depth)
  return row

          
def by_chromosome(table):
  """Iterate over bins grouped by chromosome name."""
  for chrom, subtable in table.groupby("chromosome", sort=False):
    yield chrom, table.as_dataframe(subtable)
    


def coverage_bedcov(bed_fname, bam_fname, min_mapq, procs=1):
  """Calculate coverage depth in parallel using pysam.bedcov
  PARALLEL: bed_fname splited into raw chunks (temp files)
  OUTPUT: pandas dataframe
          chr1 start  end depth norm  log2
  """
  if procs == 1:
    table = bedcov(bed_fname, bam_fname, min_mapq)
  else:
    chunks = []
  with futures.ProcessPoolExecutor(procs) as pool:
    args_iter = ((bed_chunk, bam_fname, min_mapq) for bed_chunk in to_chunks(bed_fname))
    for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
      chunks.append(table)
      rm(bed_chunk_fname)
    table = pd.concat(chunks, ignore_index=True)
    
  if 'gene' in table:
    table['gene'] = table['gene'].fillna('-')
  else:
    table['gene'] = '-'
    
  spans = table.end - table.start
  #ok_idx = (spans > 0)
  table = table.assign(depth=0, log2=NULL_LOG2_COVERAGE)
  #table.loc[ok_idx,'depth'] = (table.loc[ok_idx, 'basecount']/spans[ok_idx])  
  table['depth'] = table['basecount']/spans  ## DEPTH
  
  tot_region_size = spans.sum()
  tot_basecount = table.basecount.sum()
  ave_depth = tot_basecount/tot_region_size
  
  table['norm'] = table['depth']/ave_depth   ## DEPTH normalized to 1
  
  return table

def _bedcov(args):
  """Wrapper for parallel."""
  bed_fname = args[0]
  table = bedcov(*args)
  return bed_fname, table

def bedcov(bed_fname, bam_fname, min_mapq):
  """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.
  
  i.e. mean pileup depth across each region.
  """
  # Count bases in each region; exclude low-MAPQ reads
  cmd = [bed_fname, bam_fname]
  if min_mapq and min_mapq > 0:
    cmd.extend(['-Q', bytes(min_mapq)])
  try:
    raw = pysam.bedcov(*cmd, split_lines=False)
  except pysam.SamtoolsError as exc:
    raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
  if not raw:
    raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
  columns = detect_bedcov_columns(raw)
  table = pd.read_csv(StringIO(raw), sep='\t', names=columns, usecols=columns)
  return table

def detect_bedcov_columns(text):
    """Determine which 'bedcov' output columns to keep.

    Format is the input BED plus a final appended column with the count of
    basepairs mapped within each row's region. The input BED might have 3
    columns (regions without names), 4 (named regions), or more (arbitrary
    columns after 'gene').
    """
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount < 3:
        raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
    if tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    if tabcount == 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
    # Input BED has arbitrary columns after 'gene' -- ignore them
    fillers = ["_%d" % i for i in range(1, tabcount - 3)]
    return ['chromosome', 'start', 'end', 'gene'] + fillers + ['basecount']
