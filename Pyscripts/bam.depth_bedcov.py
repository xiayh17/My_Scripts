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


def main(bed_fname, bam_fname, by_count, min_mapq, processes):
  if by_count:
    #TODO
  elif by_bedcov:
    coverage_bedcov(bed_fname, bam_fname, by_count, min_mapq, processes)

      

    
    
    
    


def coverage_bedcov(bed_fname, bam_fname, by_count, min_mapq, processes):
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
