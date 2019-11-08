'''
pysam pileup USAGE, compared with samtools mpileup
'''

import sys
import pysam
#samfile = pysam.AlignmentFile(sys.argv[1], "rb" )
samfile = pysam.AlignmentFile(sys.argv[1], "rb" )
#for pileupcolumn in samfile.pileup(ignore_overlaps=True,flag_filter=False,ignore_orphans=True,min_base_quality=False,min_mapping_quality=False):
#for pileupcolumn in samfile.pileup(ignore_overlaps=False,flag_filter=False,ignore_orphans=True,min_base_quality=False,min_mapping_quality=False): ## samtools mpileup -Q 0 -q 0 --count-orphans
#for pileupcolumn in samfile.pileup(ignore_overlaps=False,flag_filter=False,ignore_orphans=False,min_base_quality=False,min_mapping_quality=False): ## samtools mpileup -Q 0 -q 0

#for pileupcolumn in samfile.pileup(ignore_overlaps=False,flag_filter=0,ignore_orphans=True,min_base_quality=False,min_mapping_quality=False):  ## samtools  mpileup -Q 0 -q 0 -ff 0
#for pileupcolumn in samfile.pileup(ignore_overlaps=False,flag_require=0,ignore_orphans=True,min_base_quality=False,min_mapping_quality=False): ## samtools  mpileup -Q 0 -q 0 -rf 0
#for pileupcolumn in samfile.pileup(ignore_overlaps=False,ignore_orphans=True,min_base_quality=False,min_mapping_quality=False):
    print ("coverage at base %s = %s" %
           (pileupcolumn.pos+1, pileupcolumn.n))
    #for pileupread in pileupcolumn.pileups:
        #if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            #print ('\tbase in read %s = %s' %
             #     (pileupread.alignment.query_name,
             #      pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()
