### Code from https://bioinformatics.stackexchange.com/questions/2242/how-to-convert-bed-to-gff3
#The .tsv contians chromosome names and lengths in two columns
import os
import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bed", help = "specify a bed file to use as input (bed12 format)")
parser.add_argument("-c", "--chrom", help = "specify a .tsv input file containing chromosome names and lengths in two columns")
parser.add_argument("-o", "--output", help = "specify name of the .gtf output file ")
args = parser.parse_args()

try:
   args.bed
   args.chrom
   args.output
   a = os.path.isfile(args.bed)
   b = os.path.isfile(args.chrom)
   if a and b == True :
      ### This is the code that does the actual processing of the .bed into .gtf
      print("Processing...")
      with open(args.output, "a") as f:  
         for line in open(args.chrom):
            fields = line.strip().split("\t")
            print(fields[0], ".", "contig","1",str(fields[1]), ".", "+", ".", "ID=%s" % fields[0], file = f, sep='\t')

         for line in open(args.bed):
            fields = line.strip().split("\t")
            contig = fields[0]
            # note: BED is 0-based, half-open, GFF is 1-based, closed
            start= (int(fields[1]) + 1)
            end=fields[2]
            name=fields[3]
            score=fields[4]
            strand=fields[5]
            print(contig, "from_file:" + args.bed, "element", str(start), end, score, strand, ".", "ID=%s;parent=%s" % (name, contig), sep='\t', file=f)

            block_sizes = map(int, fields[10][:-1].split(","))   # The [:-1] trims the last element from the field which is a "," If this is not trimmed an additional '' value exists in the list
            block_starts = map(int,fields[11][:-1].split(","))    # which int() cannot act on, resulting in a traceback
         
            for (block, (bstart, blen)) in enumerate(zip(block_starts, block_sizes)):
               bend = start - 1 + bstart + blen   # without the -1 the bend extends 1 position after the element's end
               print(contig, "from_file:" + args.bed, "sub-element", str(start + bstart), str(bend), score, strand, ".", "ID=%s_%i;parent=%s" %(name, block, name), sep='\t', file = f)
      print("GTF file " + args.output +" completed!")
      ###
   else:
      print("Invalid inputs; check paths and ensure parameters are given for --bed --chrom and --output")
except:
   print("Invalid inputs; check paths and ensure parameters are given for --bed --chrom and --output")
