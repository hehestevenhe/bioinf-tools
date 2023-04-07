# Python script to convert bed12 file into gtf format 

For manual input prompts, run using the below command with Python3 installed.  
  
The user will be prompted for a .tsv file containing chromosome names and sizes (in two columns), as well as the bed12 file and a desired output name.  
<pre><code>python bed12-to-gtf.py </pre></code>

Alternatively, another version of the script has also been provided which allows the arguments to be passed directly from command line using the code below

<pre><code>python bed12-to-gtf-CLI.py -b [bed file] -c [tsv file] -o [output name] </pre></code>

The following files are included here:
  - The python script with input prompts (bed12-to-gtf.py) or command line input (bed12-to-gtf-CLI.py)
  - An example .tsv file (downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)
  - An example bed12 file (10 lines)
  - An example output gtf file using the two input files above
