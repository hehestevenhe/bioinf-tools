### Code adapted from https://bioinformatics.stackexchange.com/questions/2242/how-to-convert-bed-to-gff3
### The intended bed file inputs originally come from gtf files which have been converted using bedparse::gtf2bed and processed in TRAWLING
### This gtf2bed method only parses transcript, exon and CDS/start_codon/stop_codon information from the gtf file; all other features are lost (eg. UTRs)
### As a result, the output from this script only contains transcript, exon and CDS features (start and stop codon info is unretrievable from the gtf2bed conversion)
### Code for this initial conversion which was helpful for reversing the conversion can be found at https://github.com/tleonardi/bedparse/blob/master/bedparse/converters.py
import os
import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bed", help = "specify a bed file to use as input (bed12 format)")
parser.add_argument("-o", "--output", help = "specify name of the .gtf output file ")
args = parser.parse_args()

try:
   args.bed
   args.output
   a = os.path.isfile(args.bed)
   if a == True :
      ### This is the code that does the actual processing of the .bed into .gtf
      print("Processing...")
      with open(args.output, "a") as f:  
         for line in open(args.bed):
            fields = line.strip().split("\t")
            contig = fields[0]
            # note: BED is 0-based, half-open, GTF is 1-based, closed
            # GTF CDS annotation include the start but NOT the stop codon (https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md)
            start=(int(fields[1]) + 1)
            end=fields[2]
            name=fields[3]
            score=fields[4]
            strand=fields[5]
            cds_start=(int(fields[6]) + 1)
            cds_end=(int(fields[7]) - 3)
            gene_name=fields[12]    ### These are extraFields from gtf2bed which have no strict convetion in bed 
            gene_id=fields[13]      ### make sure they are present and column numbers are correct!
            transcript_biotype=fields[14]
            print(contig, "from_file:" + args.bed, "transcript", str(start), end, score, strand, ".", 'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_biotype "%s";' % (gene_id, gene_name, name, transcript_biotype), sep='\t', file=f)
            if fields[6]!=fields[7]:
               print(contig, "from_file:" + args.bed, "CDS", str(cds_start), cds_end, score, strand, ".", 'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_biotype "%s";' % (gene_id, gene_name, name, transcript_biotype), sep='\t', file=f)
               print(contig, "from_file:" + args.bed, "start_codon", str(cds_start), (cds_start + 2), score, strand, ".", 'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_biotype "%s";' % (gene_id, gene_name, name, transcript_biotype), sep='\t', file=f)
               print(contig, "from_file:" + args.bed, "stop_codon", (cds_end + 1), (cds_end + 3), score, strand, ".", 'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_biotype "%s";' % (gene_id, gene_name, name, transcript_biotype), sep='\t', file=f)
            block_sizes = map(int, fields[10][:-1].split(","))    # The [:-1] trims the last element from the field which is a "," If this is not trimmed an additional '' value exists in the list
            block_starts = map(int,fields[11][:-1].split(","))    # which int() cannot act on, resulting in a traceback
         
            for (block, (bstart, blen)) in enumerate(zip(block_starts, block_sizes)):
               bend = start - 1 + bstart + blen   # without the -1 the bend extends 1 position after the element's end
               print(contig, "from_file:" + args.bed, "exon", str(start + bstart), str(bend), score, strand, ".", 'gene_id "%s"; gene_name "%s"; ID "%s_%i"; transcript_id "%s"; transcript_biotype "%s";' %(gene_id, gene_name, name, block, name, transcript_biotype), sep='\t', file = f)
      print("GTF file " + args.output +" completed!")
      ###
   else:
      print("Invalid inputs; check paths and ensure parameters are given for --bed and --output")
except:
   print("Invalid inputs; check paths and ensure parameters are given for --bed and --output")
