### Adjust the following input parameters as required
############################################################################
input_lncrna_ref = "reference-annotation/lncipedia_5_2_hg38.gtf"           # A reference gtf file specifying the genomic coordinates for your lncRNAs of interest; gene names must be present in the attributes column
lncrna_header = 2                                                          # The number of header/comment rows prior to gtf entries in "input_lncrna_ref"
input_coding_ref = "reference-annotation/hg38.gtf"                         # A reference gtf file containing the genomic coordinates of protein coding genes e.g. hg38 reference gtf
coding_header = 5                                                          # The number of header/comment rows prior to gtf entries in "input_coding_ref"
input_cnv_ref = "reference-annotation/CCLE_ABSOLUTE_combined_20181227.txt" # Reference data downloaded from depmap containing CNV annotation "https://depmap.org/portal/download/all/?release=CCLE+2019&file=CCLE_ABSOLUTE_combined_20181227.xlsx"                                                                      
input_genes = "./example-gene-input.tsv"                                   # A tsv file with one column named "gene" listing all lncRNA genes of interest
input_depmap_ID = "./example-depmap_ID-input.tsv"                          # A tsv file with one column named "depmap_ID" listing the depmap IDs for all cell lines of interest
output_prefix= "lncRNA"                                                    # Prefix to be specified for any output files generated
############################################################################

import pandas as pd
import re

gtf_init=[]
with open(input_lncrna_ref) as gtf:
    for _ in range(lncrna_header):
        next(gtf)
    for line in gtf:
        gtf_init.append(line)
ref_gtf=pd.DataFrame([entry.strip().split('\t') for entry in gtf_init], 
                     columns=('contig', 'source','feature','start','end','score','strand','frame','attribute'))


ref_init=[]
with open(input_coding_ref) as ref:
    for _ in range(coding_header):
        next(ref)
    for line in ref:
        temp = line.strip().split('\t')
        if temp[2] == "gene":
            ref_init.append(temp)
        else:
            pass
ref_coding=pd.DataFrame(ref_init,columns=('contig', 'source','feature','start','end','score','strand','frame','attribute'))

cnv_ref = pd.read_csv(input_cnv_ref, sep = '\t')

gene_list = pd.read_csv(input_genes,sep = '\t')
depmap_ID_list = pd.read_csv(input_depmap_ID, sep = '\t')

missing_genes = 0
total_overlap_genes = pd.DataFrame(columns=['input_gene','contig','source','feature','start','end','score','strand','frame','attribute'])
total_overlap_cnv = pd.DataFrame(columns=['input_gene','sample',
 'Chromosome',
 'Start',
 'End',
 'Num_Probes',
 'Length',
 'Modal_HSCN_1',
 'Modal_HSCN_2',
 'Modal_Total_CN',
 'Subclonal_HSCN_a1',
 'Subclonal_HSCN_a2',
 'Cancer_cell_frac_a1',
 'Ccf_ci95_low_a1',
 'Ccf_ci95_high_a1',
 'Cancer_cell_frac_a2',
 'Ccf_ci95_low_a2',
 'Ccf_ci95_high_a2',
 'LOH',
 'Homozygous_deletion',
 'depMapID'])
for x in range(len(gene_list)):
    gene = gene_list.iloc[x]['gene']
    try:
        # Generates start and end positions for the lncRNA gene
        import re
        df=ref_gtf[ref_gtf['attribute'].str.contains(gene)]
        gene_chrom = re.findall(r'\d+', df.iloc[0]['contig'])
        gene_start = int(df.iloc[0]['start'])
        gene_end = int(df.iloc[-1]['end'])
        # Checks whether lncRNA overlaps with any other genes
        overlaps_init = []
        for i in range(len(ref_coding)):
            if ref_coding['contig'][i] == gene_chrom[0]:
                if gene_start <= int(ref_coding['start'][i]) <= gene_end:
                    overlaps_init.append(ref_coding.iloc[i])
                else:
                    pass
            else:
                pass
            if ref_coding['contig'][i] == gene_chrom[0]:
                if gene_start <= int(ref_coding['end'][i]) <= gene_end:
                    overlaps_init.append(ref_coding.iloc[i])
                else:
                    pass
            else:
                pass
        overlaps = pd.DataFrame(overlaps_init)
        overlaps = pd.DataFrame.drop_duplicates(overlaps)
        overlaps.insert(0,"input_gene", gene)
        total_overlap_genes=pd.concat([total_overlap_genes,overlaps])
        
        aggregate_overlaps=pd.DataFrame(columns=list(cnv_ref.columns))
        for y in range(len(depmap_ID_list)):
            depmap_ID = depmap_ID_list.iloc[y]['depmap_ID']
            # Checks whether lncRNA overlaps with known CNV region
            cnv_init = []
            filter_cnv = cnv_ref[cnv_ref['depMapID']==depmap_ID]
            for i in range(len(filter_cnv)):
                if filter_cnv.iloc[i]['Chromosome'] == int(gene_chrom[0]):
                    if gene_start <= int(filter_cnv.iloc[i]['Start']) <= gene_end:
                        cnv_init.append(filter_cnv.iloc[i])
                    else:
                        pass
                    if gene_start <= int(filter_cnv.iloc[i]['End']) <= gene_end:
                        cnv_init.append(filter_cnv.iloc[i])
                    else:
                        pass
                else:
                    pass
            cnv_overlap = pd.DataFrame(cnv_init)
            cnv_overlap = pd.DataFrame.drop_duplicates(cnv_overlap)
            aggregate_overlaps=pd.concat([aggregate_overlaps,cnv_overlap])
        total_overlap_cnv=pd.concat([total_overlap_cnv,aggregate_overlaps])
    except:
        print("Warning: " + gene + " could not be found within the LNCipedia reference data")
        missing_genes += 1
if total_overlap_genes.empty:
    print("No overlapping genes were identified against input genes")
else:
    pd.DataFrame.to_csv(total_overlap_genes, sep='\t', path_or_buf= output_prefix+"_overlapping_genes.tsv", index=False)
print(str(missing_genes) + " input gene(s) could not be found within the LNCipedia reference data")
if total_overlap_cnv.empty:
    print("No overlapping CNV regions were identified against input genes for all input depmap IDs")
else:
    pd.DataFrame.to_csv(cnv_overlap, sep='\t', path_or_buf= output_prefix+"_overlapping_cnv.tsv", index=False)