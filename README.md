This GitHub archive provides a collection of scripts used to undertake various bioinformatic tasks described in the manuscript ‘Target enrichment of long open reading frames and ultraconserved elements to link microevolution and macroevolution in non-model organisms’ by Ortiz-Sepulveda et al. (2022). Specifically, these scripts are part of a workflow to select and jointly enrich long open reading frames (ORFs) and ultraconserved elements (UCEs) from genomic samples. The workflow is summarized in Fig. 1 of Ortiz-Sepulveda et al. (2022) and they cover the transcriptomic steps for the selection of ORFs, the comparative genomics for the selection of UCEs, up to bioinformatic analyses concerned with probe design and the sequencing of genomic libraries. Below we outline each script and in what task it has been used.

# 1 Identifying target open reading frames
## 1.1 Transcriptome sequencing
No scripts

## 1.2 Transcriptome assembly
``` 01-Assembly_Transcriptomes.py ```

This python script takes filtered transcriptomic reads as input and allows to assemble de novo a transcriptome per individual using the AGALMA v. 2.0.0 pipeline of Dunn et al. (2013, BMC Bioinformatics, 14, e330).

``` 02-transcript_qual.sh ```

This bash script allows to perform quality control on de novo assembled transcriptomes using TransRate v.1.0.3 of Smith-Unna et al. (2016, Genome Research, 26, 1134-1144).

``` 03-busco_script.sh ```

After obtaining de novo assembled transcriptomes (either raw AGALMA transcriptomes, or  transcriptomes that underwent quality control with TransRate), we have verified transcriptome completeness by benchmarking against the BUSCO Metazoa_odb9 database (see Simão et al. 2015, Bioinformatics, 31, 3210-3212; Waterhouse et al. 2017 Molecular Biology and Evolution, 35, 543-548). This bash script indicates how such benchmarking was performed.
 
## 1.3 Selection of open reading frames as target
``` 04-cdhit_clustering.sh ```

This bash script allows to cluster several de novo assembled transcriptomes (either raw AGALMA transcriptomes, or  transcriptomes that underwent quality control with TransRate) into a supertranscriptome by comparing and clustering contigs in groups of transcripts with >95% similarity. The software used for such clustering is CD-Hit-Est v. 4.8.1 (Li and Godzik, 2006, Bioinformatics, 22, 1658-1659). 

``` 05-README_transdecoder ```

This README file outlines how ORF predictions were performed with TransDecoder v. 5.5.0 (Haas & Papanicolaou, 2018, retrieved from http://transdecoder.github.io), how the resulting ORFs were again clustered and prepared for selection of target ORFs from comparisons to the BUSCO Metazoa_odb9 database and the Unioverse probe regions of Pfeiffer et al. (2019, Molecular Phylogenetics and Evolution, 137, 114-126). 

``` 06-README_selection_ORFs_BUSCO ```

This README file outlines the general procedure that was used to select ORF targets from comparison to the BUSCO Metazoa_odb9 database. It uses the UTILITY_get_subfasta.py script to substract a .fasta sequence from a .fasta file based on the name of the contig. It also includes a script to perform alignments with MUSCLE v. 3.5.1551 (Edgar, 2004, Nucleic Acids Research, 32, 1792-1797). 

``` 07-README_selection_ORFs_Unioverse ```

This README file outlines the general procedure that was used to select ORF targets from comparison to the Unioverse probe set. It uses the UTILITY_get_subfasta.py script to substract a .fasta sequence from a .fasta file based on the name of the contig. It starts by using Yass v. 1.15 (Noé and Kucherov, 2005, Nucleic Acids Research, 33, W540-W543) to map the Unioverse probe regions to our clustered ORFs.

# 2 Identifying target ultraconserved element
``` 08-search_UCE.py ```

This python script allows to examine various input genomes to select shared ultraconserved elements. It uses the PHYLUCE pipeline v. 1.6.8 (Faircloth, 2016, Bioinformatics, 32, 786-788) to find UCEs and also performs several types of post-treatment on PHYLUCE-generated data.

# 3 Identity and overlap screening
Comparative examinations of identity and overlap screening among the ORFs and UCEs that were selected as target regions were performed with Yass. These examinations follow the use of Yass as indicated in e.g. ``` 07-README_selection_ORFs_Unioverse ```

# 4 Evaluation of retained target regions
``` 09-README_simulation_ARTilluminaReads_Venustaconcha ```

This README file outlines the general procedure that was used to generate error-free ART_Illumina reads to map to the Venustaconcha genome. It makes use of ART (ART_Illumina) v. 2.5.8 (Huang et al., 2012, Bioinformatics 28, 593-594). This procedure allows to evaluate in silico how probes for our selected ORF and UCE targets map to the distant Venustaconcha genome. After this procedure target regions were submitted to Arbor Bioscientific for probe design. 

# 5 Illumina NextSeq sequencing 
No scripts, postprocessing of these sequences was done with pipelines that have been described elsewhere, i.e. for ORF targets with HybPiper v.1.3.1 (Johnson et al. 2016, Applications in Plant Sciences, 4, e1600016) and UCE targets with Phyluce v.1.6.8 (Faircloth 2016, Bioinformatics, 32, 786-788).
