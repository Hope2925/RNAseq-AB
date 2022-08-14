# RNAseq-AB
Platform for performing and comparing differential expression RNA-seq results from EBseq &amp; DEseq analysis with applicable add-ons of GFF3 data.

Expression files can be analyzed using two different statistical methods, DESeq & EBSeq. DESeq is more readily used by the scientific community but EBSeq is better known for experiments with low biological replicates. 

## Data required before use:
  - Directory called exp3 with text files from RSEM expression analysis. Guided use of RSEM to produce appropriate results can be found in the Prep folder: "TerminalCodeRNASeq24PhrB.docx."
    -  This directory should be housed with both .Rmd files.
  - Directory called Genomes with the reference genomes (formed by GFF3 and Fasta files as described in "TerminalCodeRNASeq24PhrB.docx) and blast results as a text file (as produced by "TerminalCodeRNASeq24PhrB.docx) 
 
## EBSeq:
All differential expression analysis was done on RNASeqPhrB.Rmd file using R 4.1.2. To use this file, the directory exp3 & this document must be in the same directory.

*Note: While DESeq using adjusted p-values to determine significance, EBSeq using False Discovery Rates (FDR) and PPDE (posterior probability for DE).*
    
  ### Initial Output:
  - found in directory RNASeqPhrB_files (unless changed by user in Rmd script)
  - GeneFC**.csv = Data (PPDE, FC, expression levels) for all genes. This is what would be used as the expression file uploaded to NCBI.
  - FinalMatrixDE**.csv = Same as above but ONLY for differentially expressed genes
     
   ### Analyzed Output:
   - found in directory EBSeq
   - Final<>EB.csv = FinalMatrixDE but with BLAST information included (Protein ID, alternative locus b/n A!S and ACX60 genomes) using the python script.
      - *Note: You can color-code and sort with an excel file but be sure to save as a new excel file.*
      - ALL indicates ALL genes rather than just DE genes
   - Vennattempt = Venn diagram produced from R block Visualization

## DESeq:
All differential expression analysis was done on RNASeqDEseq.Rmd file using R 4.1.2. To use this file, the directory exp3 & this document must be in the same directory.

The comprehensive report produced by DESeq can be produced as a html file by running the Visualization block of code (lines 142-146).

   ### Output:
   - found in directory exp3
   - <>_norm.csv = normalized expression values for each gene
   - <>_results = differential expression results with log2FC, p-values, and normalized expression values for each gene
   - DESeq/DESeqBlsA.csv = _results but with BLAST info included (Protein ID, alternative locus b/n A1S and ACX60 genome) using the python script Comparing results Final.py. 


## Applicable versions/parameters:
  •	Trimming
    - Trim Galore vs 0.6.6
    - Cutadapt vs 3.5
    - Quality Phred score cut off of 20
  
  •	Reads that aligned using bowtie2 
  •	R 4.1.2

