### Code from https://bioinformatics.stackexchange.com/questions/2242/how-to-convert-bed-to-gff3
#The .tsv contians chromosome names and lengths in two columns
import os
import os.path
def program():
   chrom_tsv = input("Enter chromosome length tsv file: ")
   bed12 = input("Enter bed12 file: ")
   output = input("Please name the output .gtf file: ")
   while True:
      print("Please confirm the above inputs are correct (y/n)")
      confirm = input()
      if confirm in ("y","n"):
         break
      print("Invalid input: Please enter 'y' or 'n' ")
            
   if confirm == "y":
      a = os.path.isfile(chrom_tsv)
      b = os.path.isfile(bed12)
      if a and b == True:
         print("Processing...")
         ### This is the code that does the actual processing of the .bed into .gtf
         with open(output+".gtf", "a") as f:  
            for line in open(chrom_tsv):
               fields = line.strip().split("\t")
               print(fields[0], ".", "contig","1",str(fields[1]), ".", "+", ".", "ID=%s" % fields[0], file = f)
         
            for line in open(bed12):
               fields = line.strip().split("\t")
               contig = fields[0]
               # note: BED is 0-based, half-open, GFF is 1-based, closed
               start= (int(fields[1]) + 1)
               end=fields[2]
               name=fields[3]
               score=fields[4]
               strand=fields[5]
               print(contig, "from_file:" + bed12, "element", str(start), end, score, strand, ".", "ID=%s;parent=%s" % (name, contig), sep='\t', file=f)

               block_sizes = map(int, fields[10][:-1].split(","))   # The [:-1] trims the last element from the field which is a "," If this is not trimmed an additional '' value exists in the list
               block_starts = map(int,fields[11][:-1].split(","))    # which int() cannot act on, resulting in a traceback
            
               for (block, (bstart, blen)) in enumerate(zip(block_starts, block_sizes)):
                  bend = start - 1 + bstart + blen   # without the -1 the bend extends 1 position after the element's end
                  print(contig, "from_file:" + bed12, "sub-element", str(start + bstart), str(bend), score, strand, ".", "ID=%s_%i;parent=%s" %(name, block, name), sep='\t', file = f)
         print("GTF file " + output + ".gtf completed!")
         ###
      else:
         print("Input files cannot be found -- include file paths if these are not in the same directory as this .py file\nExiting program")
         exit()
   else:
      program()
program()
