# PlanTAMI
###### Tool for PLants Transcriptomes Analyisi with Mutual Information

This **tool** takes two lists of DEGs of different species and returns the the directo orthologs and orthologs family enriched
between them, to eavoulate the family Normalise Pointwise Mutual Information is used and an empirical p value is computed.
In this repository are available also the results obtained compraing the inflorescence transcriptome of Arabidopsis, Rice, Tomato and Barley. Also a detailed analysis of the more promisong family is provided. For detailed description of the tool and the dataset analysed refer to the file: Prototyping, implementation and validation of PlanTaMI: a tool for Plant Transcriptomics using Mutual Information. 

The **input** needed by the tool are:
  + --plaza 
  table with geneID, speciesID and familyID semicolon separed (https://bioinformatics.psb.ugent.be/plaza/)
  + --sp1 , --sp2 
  the species of DEGs (es. ath, sly, see plaza web site)
  + --in_sp1 --in_sp2
  DEGs lists with the same nomenclature of the plaza table
  
 The **optional arguments** are:
  + --th_sc
  the threshold for calling significant direct orthologs in the hypergeometric test (def=0.05)
  + --th
  the threshold for calling significant orthologs family corrected p-value in the NPMI analysis (def=0.05)
  + --random
  the number of list of rabndomID to compute empirical p-value (def=10000)
  + --FDR
  the method to adjust the p-value, Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY) are available (def=BY)
  + --sample
  the name of the sample to save results (def=my_sample_result)
  + --verbose
  the style of the output printed on the terminal during the run (0=none, 1=compact, 2=full, 3=line, def=1)

The **output** return by the tool are: 
 + -significant_direct_orthologs (or NOT_significant_direct_orthologs)
 list of direct orthologs between the two species, the fold enrichment and the relative p-value 

 + -significant_genes_sp1
 list of gene IDs of sp1 resulted significant in the NPMI analysis 
 
 + -significant_genes_sp2
 list of gene IDs od sp2 resulted significant in the NPMI analysis 
 
 + -significant_common_families_genes
 list of families significant in both sp1 and sp1 and the relative genes
 
 + -all_common_families_genes
 list of all families in common between sp1 and sp1 and the relative genes
 
 + -NPMI_results_sp1
 table with NPMI value for each family of sp1 in common with sp2, in particulare are reported:
   + Fam_ID: Family ID
   + DEGs: DEGs in the family and family size
   + DEGs_set: DEGs in family of the same size in the input set
   + DEGs_genome: DEGs in family of the same size in the genome
   + p_my: DEGs in the set/DEGs in family of the same size in the set
   + p_all: (family size * DEGs in family of the same size in the set) / DEGs in family of the same size in the genome
   + npmi: log2(p_my/p_all)/-log2(p_all)
   + p_value: number of times i see higer or equal NPMI value in random set / number of set generated 
   + FDR_bonf : Bonferroni multiple test correction 
   + FDR_BY: Benjamini-Yekutieli (or Benjamini-Hochberg) multiple test correction 

 + -NPMI_results_sp2
 as above but for sp2
 
 + -stat
 feature of the backgroud, the sample and some easy to read results 
 
  
