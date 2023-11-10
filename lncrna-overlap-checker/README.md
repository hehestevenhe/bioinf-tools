_A python-based tool to check whether lncRNA genes of interest overlap with other coding genes and/or regions of known CNV in specified cell lines_

### Reference file download
__lncRNA reference gtf file__: https://lncipedia.org/download <br>
__hg38 reference gtf file__: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ <br>
__depmap CNV annotation info__: https://depmap.org/portal/download/all/?release=CCLE+2019&file=CCLE_ABSOLUTE_combined_20181227.xlsx <br>

_Please note that the depmap CNV annotation info must be converted from .xlsx to .tsv before being used as input_

### Running the code
Input variables are denoted in the code and must be manually edited, with comments providing additional details. The code has been provided in both .ipynb and .py formats, with the Jupyter notebook providing more information on each code module. <br>

The code generates 2 potential .tsv files as outputs for any given gene:
1.  _"overlapping_genes.tsv"_ returns the gtf entries of genes that overlap with the genomic coordinates of the lncRNA gene of interest. Note: A new column specifies which gtf entries pertain to any given input gene.

2.  _"overlapping_cnv.tsv"_ return the entries of CNV regions that overlap with the genomic coordinates of the lncRNA gene of interest, keeping with the format of the original depmap CNV reference annotation. Note: A new column specifies which entries pertain to any given input gene.

__These outputs are not generated if no relevant overlaps are detected__ 

In the case that the lncRNA gene of interest cannot be found in the lncRNA reference data, a warning message will appear indicating the specific lncRNA gene

