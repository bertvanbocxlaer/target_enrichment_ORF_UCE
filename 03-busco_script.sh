#!/bin/bash

module load busco

tr=$(ls ../good_assemblies/*.fa)

for f in $tr
do 
	out=$(echo $f | cut -d '.' -f3 | cut -d '/' -f3)
	cmd="run_busco -i $f -c 20 -o BUSCO_$out -l /eep/projects1/Evolink/Reference_transcriptomes/BUSCO_DB/metazoa_odb9 -m tran"
	$cmd
done

