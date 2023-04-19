# Python script to convert bed12 file into gtf format 
The bed12-to-gtf.py python script converts a bed12 format file to gtf, with the intended bed file having been generated using TRAWLING (https://github.com/CCI-CDDB/trawling)    

From the input bed, it will generate a gtf file with the following characterisitics:  
  - Features (column 3): trascript, exon, CDS, start_codon and stop_codon features only, based on the bed input
  - Attributes (column 9): gene_id, gene_name, transcript_id and transcript_biotype, provided these are included in the bed file

**Please note:** The script assumes that gene_name, gene_id and transcript_biotype are formatted respectively in the 13th, 14th and 15th columns of the bed file  

This script takes two arguments (an input .bed file and desired output name) which are passed directly from command line following the usage below

<pre><code>python bed12-to-gtf.py -b [bed file] -o [output name] </pre></code>

The following files are included here:
  - The python script (bed12-to-gtf.py)
  - An example bed12 file (test.bed)
  - An example output gtf file using the test input above (outputtest.gtf)
