
#!/bin/bash

#module load cdhit/4.8.1

fasta=$(ls *.fa)
#outdir=$clustered_good_assemblies

for f in $fasta
do
	out=$(echo $f | cut -d '.' -f1)
	cd-hit-est -c 0.95 -M 8000 -T 15 -i $f -o clustered_good_assemblies/$out-good_clustered_95p.fa
done
