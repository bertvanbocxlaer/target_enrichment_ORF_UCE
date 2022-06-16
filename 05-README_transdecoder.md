
### Run transdecoder on either the AGALMA clustered transcriptome assemblies or the TransRate clustered assemblies
### example for TransRate clustered transcriptome assemblies

cp pipeline_Assembly/All_fasta_Transcriptomes/good_assemblies/clustered_good_assemblies/All_Unionidae_transrate_clustered_95p.fa Reference_transcriptomes/Unionidae/All_Unionidae/All_Unionidae_transrate_clustered_95p.fa
cp pipeline_Assembly/All_fasta_Transcriptomes/good_assemblies/clustered_good_assemblies/All_Unionidae_transrate_clustered_95p.fa.clstr Reference_transcriptomes/Unionidae/All_Unionidae/All_Unionidae_transrate_clustered_95p.fa.clstr

module load transdecoder

TransDecoder.LongOrfs -t All_Unionidae_transrate_clustered_95p.fa --output_dir AllUnionidae_TransRate_BestORF
TransDecoder.Predict -t All_Unionidae_transrate_clustered_95p.fa --output_dir AllUnionidae_TransRate_BestORF --single_best_only

cp All_Unionidae_transrate_clustered_95p.fa.transdecoder.cds All_Unionidae_transrate_ORFs.fa

### For TransRate clustered transcriptome assemblies we obtained 133,395 ORFs starting from 701,510 contigs
### For AGALMA transcriptome assemblies we obtained 245,861 ORFs starting from 988,460 contigs

module load busco

run_busco -i All_Unionidae_transrate_ORFs.fa -c 20 -o BUSCO_Unionidae_TransRate_ORFs -l /eep/projects1/Evolink/Reference_transcriptomes/BUSCO_DB/metazoa_odb9 -m tran

module load cdhit

nohup cd-hit-est -c 0.95 -M 8000 -T 15 -i All_Unionidae_transrate_ORFs.fa -o All_Unionidae_transrate_ORFs_95p.fa

### at this stage we have 99,359 ORFs that remain for the TransRate filtered assemblies and 131,503 ORFs for the AGALMA transcriptome assemblies

run_busco -i All_Unionidae_transrate_ORFs_95p.fa -c 20 -o BUSCO_Unionidae_TransRate_ORFs_95p -l /eep/projects1/Evolink/Reference_transcriptomes/BUSCO_DB/metazoa_odb9 -m tran

