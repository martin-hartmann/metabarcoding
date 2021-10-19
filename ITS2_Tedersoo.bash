# to start the pipeline, move into the project directory and create subdirectories "1_raw", "2_processed"
## 1_raw contains the raw fastq file (gz compressed)

# PREPARE
cd 1_raw
parallel 'gunzip -k {}' ::: *gz # unpack
xargs -a rename.txt -n 2 mv # rename files based on 2-column find-replace list saved in rename.txt (optional)
parallel "awk '{if (NR%4==1) {sub(\"_R1.fastq\",\"\",FILENAME); sub(\"_R2.fastq\",\"\",FILENAME); print \$1 \";sample=\" FILENAME \";\"} else {print \$0}}' {} > {}.tmp; mv {}.tmp {}" ::: *fastq # add filename to sequence header
cat *R1.fastq > ../2_processed/1.all.raw.R1.fastq # concatenate R1
cat *R2.fastq > ../2_processed/1.all.raw.R2.fastq # concatenate R2
rm *fastq # delete single files
cd ../2_processed # move to new directory

# QC
fastqc 1.all.raw.R1.fastq # QC on R1
fastqc 1.all.raw.R2.fastq # QC on R2
vsearch --fastq_eestats 1.all.raw.R1.fastq --output 1.all.raw.R1.eestats.txt # QC stats R1
vsearch --fastq_eestats 1.all.raw.R2.fastq --output 1.all.raw.R2.eestats.txt # QC stats R2

# QUALITY FILTER AND ASV/OTU DELINEATION
## databases required for phix removal and taxonmic classification deposited in "/databases/"
### phiX genome can be retrieved from https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1
### the UNITE database can be downloaded from https://unite.ut.ee/repository.php
bowtie2 -x /databases/phix -1 1.all.raw.R1.fastq -2 1.all.raw.R2.fastq --un-conc all.nophix --fast > /dev/null # remove phix
mv all.1.nophix 2.all.nophix.R1.fastq & mv all.2.nophix 2.all.nophix.R2.fastq # rename
cutadapt -g ^CANCGATGAAGAACGYRG -G ^CCTSCSCTTANTDATATGC -e 0.12 -m 1 --trimmed-only -o 3.all.trim.R1.fastq -p 3.all.trim.R2.fastq 2.all.nophix.R1.fastq 2.all.nophix.R2.fastq # trim primers
vsearch --fastq_mergepairs 3.all.trim.R1.fastq --reverse 3.all.trim.R2.fastq --fastqout 4.all.merge.notrim.fastq --fastaout 4.all.merge.notrim.fasta --fastq_truncqual 7 --fastq_allowmergestagger --fastq_minovlen 30 # merge PE reads (no length trimming to assess primer dimer quantities)
vsearch --fastq_mergepairs 3.all.trim.R1.fastq --reverse 3.all.trim.R2.fastq --fastqout 4.all.merge.fastq --fastaout 4.all.merge.fasta --fastq_truncqual 7 --fastq_allowmergestagger --fastq_minovlen 30 --fastq_minmergelen 150 # merge PE reads (with target specific length trimming)
vsearch --fastq_eestats 4.all.merge.notrim.fastq --output 4.all.merge.notrim.eestats.txt # quality stats for merged reads (without trimming)
vsearch --fastq_eestats 4.all.merge.fastq --output 4.all.merge.eestats.txt # quality stats for merged reads (with trimming)
vsearch --fastq_filter 4.all.merge.fastq --fastaout 5.all.eefilter.fasta --fastq_maxee 1 # filter by maximum expected error
vsearch --derep_fulllength 5.all.eefilter.fasta --sizeout --relabel Uniq --output 6.all.uniq.fasta # dereplicate sequences
vsearch --cluster_unoise 6.all.uniq.fasta --centroids 6.all.ASV.fasta --uc 6.all.ASV.uc --relabel ASV --sizeorder --sizein --sizeout --minsize 8 # delineate into ASVs
vsearch --uchime3_denovo 6.all.ASV.fasta --nonchimeras 7.all.ASV_nochim.fasta --uchimeout 7.all.ASV_nochim.txt --uchimealns 7.all.ASV_chimaln.txt --abskew 8  --sizein --sizeout # identify and remove chimeras
ITSx -i 7.all.ASV_nochim.fasta -o 8.all.ASV_ITSx --complement F --preserve T --save_regions ITS2 --allow_single_domain F -N 2 --require_anchor HMM # verify target
seqkit common -w 80 -n 7.all.ASV_nochim.fasta 8.all.ASV_ITSx.ITS2.fasta -o 8.all.ASV_ITSx.fasta # get common sequences
seqkit sort -N -w 80 8.all.ASV_ITSx.fasta -o 8.all.ASV_ITSx.fasta # re-sort sequences by name
vsearch --cluster_size 8.all.ASV_ITSx.fasta --centroids 8.all.OTU_ITSx.fasta --uc 8.all.OTU_ITSx.uc --id 0.97 --relabel OTU --sizeorder --sizein --sizeout # cluster into OTUs
awk -i inplace -F ";size" '{print $1}' 8.all.ASV_ITSx.ITS2.fasta # reformat metaxa file
awk -i inplace -F ";size" '{print $1}' 8.all.ASV_ITSx.fasta # reformat metaxa file
awk -i inplace -F ";size" '{print $1}' 8.all.OTU_ITSx.fasta # reformat metaxa file
ITSx -i 8.all.OTU_ITSx.fasta -o 8.all.OTU_ITSx --complement F --preserve T --save_regions ITS2 --allow_single_domain F -N 2 --require_anchor HMM # verify target
vsearch --usearch_global 4.all.merge.fasta --db 8.all.ASV_ITSx.fasta --id 0.97 --maxhits 1 --maxaccepts 0 --uc 9.all.ASV_map.uc --matched 9.all.ASV_map.fasta --otutabout 9.all.ASV_map.txt # map reads
vsearch --usearch_global 4.all.merge.fasta --db 8.all.OTU_ITSx.fasta --id 0.97 --maxhits 1 --maxaccepts 0 --uc 9.all.OTU_map.uc --matched 9.all.OTU_map.fasta --otutabout 9.all.OTU_map.txt # map reads
sort -n -k1.4 -o 9.all.ASV_map.txt 9.all.ASV_map.txt # sort ASVs
sort -n -k1.4 -o 9.all.OTU_map.txt 9.all.OTU_map.txt # sort OTUs

# TAXONOMIC CLASSIFICATION
vsearch --sintax 8.all.ASV_ITSx.ITS2.fasta --db /databases/unite.v80.ITS2_HMM.sintax --tabbedout 9.all.ASV_tax.sintax.txt --sintax_cutoff 0.8 --strand plus # assign taxonomy
vsearch --sintax 8.all.OTU_ITSx.ITS2.fasta --db /databases/unite.v80.ITS2_HMM.sintax --tabbedout 9.all.OTU_tax.sintax.txt --sintax_cutoff 0.8 --strand plus # assign taxonomy
cat 9.all.ASV_tax.sintax.txt | sort -n -k1.4 | awk '{print $1 "\t" $4}' | sed 's/,/\t/g' | sed '1 i\ASV\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > 9.all.ASV_tax.sintax.rf.txt # reformat
cat 9.all.OTU_tax.sintax.txt | sort -n -k1.4 | awk '{print $1 "\t" $4}' | sed 's/,/\t/g' | sed '1 i\OTU\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > 9.all.OTU_tax.sintax.rf.txt # reformat
vsearch --usearch_global 8.all.ASV_ITSx.ITS2.fasta --db /databases/unite.v80.ITS2_HMM.sintax --id 0.97 --maxaccepts 0 --maxrejects 0 --uc_allhits --uc 9.all.ASV_tax.lca.uc --lcaout 9.all.ASV_tax.lca.txt # assign taxonomy with LCA
vsearch --usearch_global 8.all.OTU_ITSx.ITS2.fasta --db /databases/unite.v80.ITS2_HMM.sintax --id 0.97 --maxaccepts 0 --maxrejects 0 --uc_allhits --uc 9.all.OTU_tax.lca.uc --lcaout 9.all.OTU_tax.lca.txt # assign taxonomy with LCA
cat 9.all.ASV_tax.lca.txt | sort -n -k1.4 | sed 's/,/\t/g' | sed '1 i\ASV\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > 9.all.ASV_tax.lca.rf.txt # reformat taxonomy
cat 9.all.OTU_tax.lca.txt | sort -n -k1.4 | sed 's/,/\t/g' | sed '1 i\OTU\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > 9.all.OTU_tax.lca.rf.txt # reformat taxonomy
