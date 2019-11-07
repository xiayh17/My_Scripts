# My_Scripts
## Highly used scripts for data processing or plotting
* [pysam usage docs](https://kevinzjy.github.io/2018/07/27/180727-Python-pysam/)

## Careful about depth of coverage of Overlapping Paired-End Reads
* [samtools mpileup](https://www.biostars.org/p/87299/#284421) count 1x for overlaping position by default, but can count 2x using "-x" (sets overlaping bases quality to zero), so using -Q 1 to filter
* [samtools depth](https://www.biostars.org/p/323047/) counts each mate separately, means 2x, but can using sambamba depth with "--fix-mate-overlaps" to avoid double-counting
## to calculate fragment length
* For [pysam](https://www.biostars.org/p/312885/), using "record.template_length"
* For samtools, using records in the 9th column of a bam
