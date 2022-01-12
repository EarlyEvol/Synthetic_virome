#!/bin/bash
#SBATCH --job-name=virome.gamma.5.4
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100GB
#SBATCH -e slurm-%j.%x
#SBATCH -o slurm-%j.%x
#SBATCH --partition=tb

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES

source ~/.bashrc
echo " Loading conda environment"
conda activate gen2
module load ncbi-blast
module load diamond

run_name=gamma.5
threads=18

######################################
############### TODO #################
######################################
#filter blast identified viral contigs by length. (many are sub 500bp)



#######stuff to add  for virus sequence dope in########
virus_acc=''
virus_fasta=Actinobacteriophages-All.names.fasta
number_vir_seqs=5
vir_target_seq=random_virus_$number_vir_seqs.fasta
#just to deal with some seqs being in bacteria
add_seqs=10
#for excluding seqs similar to background seqs
qident=97.5
perc_cov=0.80
##################################

#######  Database variable  #########
#CAT database (and others later probably)
database=~/scratch/database/
RefSoil_names=$database/RefSoil_plasmids/DATABASE_plasmids/ref_soils_with_species.txt
CAT_db=$database/CAT_prepare_20210107/2021-01-07_CAT_database/
CATtax_db=$database/CAT_prepare_20210107/2021-01-07_taxonomy/
###################################


#######  Background stuff #########
number_seqs=10
accession_list=random_$number_seqs
ref_name=full_background.fasta
#reads_fq=bkgnd_${number_seqs}_vir_${number_vir_seqs}
#################################


####  ART perameters  ########
fold_cov=10
read_len=150
frag_len=400
frag_stder=20
seq_tech=MSv3
bknd_cov_fractions='1 .5'
###art VIRUS peramerters#####
vir_fold_cov=20
vir_cov_fractions='1 .5 .0625'
#################


####  virus scaffold annotation stuff  ####
vir_target_scaffolds=vir_target_contigs.fa
DVF_out="scaffolds.fasta_gt3000bp_dvfpred.out_file"
CAT_out="out.CAT.contig2classification.tax.out_file"
VS_out="VirSorter.tsv"
vir_ident_summary="vir_ident_summary"

### general outputs #####
read_align_stats=read_align_stats.txt


#################################
#### Define Full Pipe to Run ####
#################################
# comment out modules to skip
# modules with previously completed
#  will be automatically skipped
FULL_PIPE () {
echo "
##################################
#### Running Synthetic Virome ####
##################################
"
Build_File_Structure
Background_reads_simulation
Viral_reads_simulation
SubSample_Simulated_Reads
Assemble_Simulated_Reads
Virus_Contig_Identification
Target_Vir_Contig_Ident
Ident_Target_Vir_Reads
Compare_Vir_Ident_to_Gold_Standard
Pull_Vir_Identified_Contigs
Quantify_Aligned_Reads
Calculate_Vir_Contig_Stats
Wrap_Up
}



Build_File_Structure () {
#######  Dir setup  ######
wd=$HOME/scratch/virome/$run_name
refs_dir=refs
reads_dir=reads
asm_dir=asm
vir_dir=virus
single_vir=single_vir
annot_dir=annot
contigs_dir=contigs
mkdir -p $wd
cd $wd
mkdir $refs_dir
mkdir $reads_dir
mkdir $contigs_dir
mkdir $vir_dir
mkdir $vir_dir/$single_vir
#mkdir $annot_dir
#######################
out_file=$wd/../$run_name.out
echo $(date) > $out_file
}
##############################


Background_reads_simulation () {
echo "
############################################
#######  Simulate background reads #########
############################################
"
### add fungi and other viruses later
#  change read simulation to camisim, metasim, others?
#  to be able to control abundance
#  simmulate  metagenomic reads ART, camisim, metasim, others?
FILE="${wd}/${reads_dir}/background_reads.gzip.complete"
if [ -f "$FILE" ] ; then
   echo "
   $func already complete ($FILE present)
   moving on!
   "
else
   echo "
   Start ${func}
   $(date)"
   cd $wd/$refs_dir
   #go from a list a accessions (like RefSoil) to a set of fastas
   cat $RefSoil_names | awk '{print $1}' | sort -R | head -n $number_seqs > $accession_list
   for I in $(cat $accession_list) ; do esearch -db nucleotide -query "$I" | efetch -format fasta > $I.bact.fasta ; done

   ##add code to fetch fungal seqs ##
   ##add code to fetch virus seqs ##
   cat *.bact.fasta  *.fung.fasta *.vir.fasta > $ref_name


   cd $wd/$reads_dir 
   ~/apps/art_bin_MountRainier/art_illumina \
   -sam -i $wd/$refs_dir/$ref_name \
   -ss $seq_tech \
   -l $read_len \
   -f $fold_cov \
   -m $frag_len \
   -s $frag_stder \
   -o ./background.art.

   echo "
   Trimming Simulated background reads
   "
   java -jar ~/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
   background.art.1.fq \
   background.art.2.fq \
   background.art.trim.1.fq.gz background.art.trimUp.1.fq.gz \
   background.art.trim.2.fq.gz background.art.trimUp.2.fq.gz \
   ILLUMINACLIP:/home/earlm/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36 &&\
   touch $FILE
fi
}



Viral_reads_simulation () {
echo "
########################################################
##########  Simulate the target Virus data   ###########
########################################################
"
#simulate virome and add to metagenome reads
##################
### add module ###
#if provided NCBI Accessions
#for I in virus_acc do;
#esearch -db nucleotide -query "$I" | efetch -format fasta > $I.fasta ; done
#cat *.fasta  $vir_target_seq

#if provided multifasta
#pull $number_vir_seqs + $add_seqs to buffer against seqs found in background
FILE="${wd}/${reads_dir}/virus_reads.gzip.complete"
if [ -f "$FILE" ] ; then
   echo "
   $func already complete ($FILE present)
   moving on!
   "
else
   echo "

   ############################################################
   Start ${func}
   $(date)
   ############################################################"
   cd $wd/$vir_dir
   cat $wd/../$virus_fasta | \
   grep ">" | sed 's/>//g' | \
   awk '{print $1}' | \
   sort -R | \
   head -n $(expr $number_vir_seqs + $add_seqs) \
   > $wd/$vir_dir/vir_accession_list.buffer

   ~/apps/bbmap/filterbyname.sh include=T \
   names=vir_accession_list.buffer \
   in=$wd/../$virus_fasta \
   out=$wd/$vir_dir/${vir_target_seq%.*}.buffer.fasta


   ##### make sure there are no close matches in background
   #align vir_seqs to background
   cd $wd/$refs_dir

   makeblastdb -in $wd/$refs_dir/$ref_name -dbtype nucl
   blastn -db $wd/$refs_dir/$ref_name \
   -query $wd/$vir_dir/${vir_target_seq%.*}.buffer.fasta \
   -num_threads $threads \
   -outfmt '6 qaccver saccver pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore' \
   > $wd/$vir_dir/$vir_target_seq.buffer.check.blst6


   cd $wd/$vir_dir
   #filter blast hits
   ######NEED TO ADD pident and perc_cov to variable to be passed to awk...has wierd syntax....
   cat $wd/$vir_dir/$vir_target_seq.buffer.check.blst6 | \
   awk '{if ($3 >= 90 && $4/$5 >= 0.50)  print $0}' \
   > $wd/$vir_dir/$vir_target_seq.buffer.check.blst6.filter
   cat $wd/$vir_dir/$vir_target_seq.buffer.check.blst6.filter | awk '{print $1}' \
   > $wd/$vir_dir/$vir_target_seq.buffer.check.blst6.filter.names

   #extract viral seq names with no blast hits in background
   cat $wd/$vir_dir/vir_accession_list.buffer | \
   grep -vf $wd/$vir_dir/$vir_target_seq.buffer.check.blst6.filter.names | \
   sort -r | head -n $number_vir_seqs \
   > vir_accession_list.clean

   #check to make sure there are enough viral seqs not found in background...
   #this is stupid and messy.
   number=$(cat vir_accession_list.clean | wc -l)
   if ! (($number >= $number_vir_seqs)) ; then
      echo "NOT ENOUGH UNIQ VIRAL SEQS, RERUN"
      exit
   fi

   #make clean viral multifasta
   ~/apps/bbmap/filterbyname.sh include=T \
   names=vir_accession_list.clean \
   in=$wd/../$virus_fasta \
   out=$wd/$vir_dir/$vir_target_seq

   #split multifasta into ingle files for running through metaquast later
   cd $wd/$vir_dir/$single_vir/
   cat $wd/$vir_dir/$vir_target_seq |\
   awk '{
      if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
      print $0 > filename
      }'

   cd $wd/$reads_dir
   #simmulate  metaVIROME reads ART, camisim, metasim, others?
   #use same parameters as background except coverage
   cd $wd/$reads_dir/
   ~/apps/art_bin_MountRainier/art_illumina \
   -sam -i $wd/$vir_dir/$vir_target_seq \
   -ss $seq_tech \
   -l $read_len \
   -f $vir_fold_cov \
   -m $frag_len \
   -s $frag_stder \
   -o ./target_virus.art.
   mv target_virus.art..sam target_virus.art.sam

   echo	"
   Trimming Simulated viral reads
   "
   java -jar ~/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
   target_virus.art.1.fq \
   target_virus.art.2.fq \
   target_virus.art.trim.1.fq.gz target_virus.art.trimUp.1.fq.gz \
   target_virus.art.trim.2.fq.gz target_virus.art.trimUp.2.fq.gz \
   ILLUMINACLIP:/home/earlm/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36 &&\
   touch $FILE
fi

}




SubSample_Simulated_Reads () {
echo "
###################################################
############ SubSample Simulated Reads ############
###################################################
"
#I could make this more efficient by only subsampling once per fraction if there are multiple fractions for both background and virus
cd $wd
for I in $bknd_cov_fractions;
do
   for J in $vir_cov_fractions;
   do
      func="Subsample reads. Fractions: $I of background and $J of viral "
      FILE="${wd}/subsmpl-bknd-${I}-vir-${J}/$reads_dir/bknd-${I}-vir-${J}.complete"
      if [ -f "$FILE" ] ; then
         echo "
         $func already complete ($FILE present)
         moving on!
         "
      else
         echo "
         Start ${func}
         $(date)"
         cd $wd
         mkdir subsmpl-bknd-${I}-vir-${J}
         cd subsmpl-bknd-${I}-vir-${J}
         mkdir $reads_dir
         cd $reads_dir

         #subsample background
         reformat.sh in1=$wd/$reads_dir/background.art.trim.1.fq.gz in2=$wd/$reads_dir/background.art.trim.2.fq.gz \
         out1=subsmplbkgnd-${I}-1.fq.tmp out2=subsmplbkgnd-${I}-2.fq.tmp samplerate=$I &&\
         mv subsmplbkgnd-${I}-1.fq.tmp subsmplbkgnd-${I}-1.fq &&\
         mv subsmplbkgnd-${I}-2.fq.tmp subsmplbkgnd-${I}-2.fq

         reformat.sh in=$wd/$reads_dir/background.art.trimUp.1.fq.gz \
         out=subsmplbkgnd-${I}-Up-1.fq samplerate=$I
         reformat.sh in=$wd/$reads_dir/background.art.trimUp.2.fq.gz \
         out=subsmplbkgnd-${I}-Up-2.fq samplerate=$I

         #subsample virus
         reformat.sh in1=$wd/$reads_dir/target_virus.art.trim.1.fq.gz in2=$wd/$reads_dir/target_virus.art.trim.2.fq.gz \
         out1=subsmplvir-${J}-1.fq.tmp out2=subsmplvir-${J}-2.fq.tmp samplerate=$J &&\
         mv subsmplvir-${J}-1.fq.tmp subsmplvir-${J}-1.fq &&\
         mv subsmplvir-${J}-2.fq.tmp subsmplvir-${J}-2.fq
         reformat.sh in=$wd/$reads_dir/target_virus.art.trimUp.1.fq.gz \
         out=subsmplvir-${J}-Up-1.fq samplerate=$I
         reformat.sh in=$wd/$reads_dir/target_virus.art.trimUp.2.fq.gz \
         out=subsmplvir-${J}-Up-2.fq samplerate=$I

         cat subsmplbkgnd-${I}-1.fq subsmplvir-${J}-1.fq \
         > bknd-${I}-vir-${J}-1.fq.tmp &&\
         mv bknd-${I}-vir-${J}-1.fq.tmp bknd-${I}-vir-${J}-1.fq

         cat subsmplbkgnd-${I}-2.fq subsmplvir-${J}-2.fq \
         > bknd-${I}-vir-${J}-2.fq.tmp &&\
         mv bknd-${I}-vir-${J}-2.fq.tmp bknd-${I}-vir-${J}-2.fq

         cat subsmplbkgnd-${I}-Up-1.fq subsmplvir-${J}-Up-1.fq \
         subsmplbkgnd-${I}-Up-2.fq subsmplvir-${J}-Up-2.fq \
         > Up-bknd-${I}-vir-${J}-0.fq.tmp && \
         mv Up-bknd-${I}-vir-${J}-0.fq.tmp Up-bknd-${I}-vir-${J}-0.fq

         if [ -f bknd-${I}-vir-${J}-1.fq ] && [ -f bknd-${I}-vir-${J}-2.fq ] && [ -f Up-bknd-${I}-vir-${J}-0.fq ];
         then
            touch $FILE
         fi
      fi
   done
done

}



Assemble_Simulated_Reads () {
echo "
###################################################
############  Assemble Simulated Reads ############
###################################################
"
cd $wd
for I in subsmpl-*;
do
   cd $wd/$I
   # assembly probably metaspades is best, could add multiple
   func="Simulated_read_assembly of $I"
   FILE="${wd}/${I}/$asm_dir/assembly.complete"
   if [ -f "$FILE" ] ; then
      echo "
      $func already complete ($FILE present)
      moving on!
      "
   else
      if [ -s ../$reads_dir/Up-bknd-*-vir-*-0.fq ];then
          single="-s ../$reads_dir/Up-bknd-*-vir-*-0.fq"
      else
          single=''
      fi
      echo "
      Start ${func}
      $(date)"
      cd $wd/$I
      mkdir asm
      cd asm
      ~/apps/SPAdes-3.15.2-Linux/bin/metaspades.py \
      -t $threads \
      -1 ../$reads_dir/bknd-*-vir-*-1.fq \
      -2 ../$reads_dir/bknd-*-vir-*-2.fq \
      $single \
      -o ./ > spades_log.txt 2>&1 &&\
      touch $FILE
   fi
done

}



Virus_Contig_Identification () {
echo "
########################################
####   virus contig identification   ###
########################################
"
cd $wd
for I in subsmpl-*;
do
   echo "######################################
   Doing viral contig identification on $I assemblies"
   cd $wd/$I
   if [ ! -d $annot_dir ]; then
	 mkdir $annot_dir
   fi
   ##################################################################################
   #Use fancy pants deeplearning program to find viral sequences in assembles contigs
   #deepvirfinder
   func="DeepVirFinder_run of $I"
   #in args
   ASM=$asm_dir/scaffolds.fasta
   DVF_threads=`expr $threads - 10`
   DVF_outdir=$annot_dir/
   #out files
   FILE=${wd}/${I}/${annot_dir}/$DVF_out
   echo "#################################"
   if [ -f "$FILE" ]; then
      echo "
      $func already complete ($FILE present)
      moving on!
      "
   else
      echo "
      Start ${func}
      $(date)"
      cd $wd/$I
      python ~/apps/DeepVirFinder/dvf.edit.py \
      -i $ASM \
      -o $DVF_outdir -l 3000 -c $DVF_threads &&\
      mv $annot_dir/scaffolds.fasta_gt3000bp_dvfpred.txt \
      $FILE
   fi
   ##################################################


   ###################################
   #run virsorter 
   #  Add checkV at some point (Already installed)
   func="Classify contigs with VirSorter $I"
   #in args
   ASM=$asm_dir/scaffolds.fasta
   VS_threads=$threads
   VS_outdir=$annot_dir/VS
   #out files
   FILE=${wd}/${I}/${annot_dir}/$VS_out
   FILE_fasta=$FILE.fasta
   echo "#################################"
   if [ -f "$FILE" ]; then
      echo "
      $func already complete ($FILE present)
      moving on!
      "
   else
      echo "
      Start ${func}
      $(date)"
      cd $wd/$I
      if [ -d annot/VS/ ]; then
         mv $VS_outdir ${VS_outdir}.bk
      fi
      virsorter run -w  $VS_outdir -i $ASM -j $VS_threads && \
      cp $VS_outdir/final-viral-boundary.tsv $FILE && \
      cp $VS_outdir/final-viral-combined.fa $FILE_fasta && \
      rm annot/VS/*
   fi
   ######################################################


   #######################################################
   #Use a huge database search to identify viral sequences
   #CAT
   #in args
   func="CAT_contig_claccification $I"
   ASM=$asm_dir/scaffolds.fasta
   CAT_threads=$threads
   CAT_db=$CAT_db
   CATtax_db=$CATtax_db
   #CAT_outdir=$annot_dir/CAT/
   #out file
   FILE=${wd}/${I}/${annot_dir}/$CAT_out
   echo "#################################"
   if [ -f "$FILE" ]; then
      echo "
      $func already complete ($FILE present)
      moving on!
      "
   else
      echo "
      Start ${func}
      $(date)"
      cd $wd/$I/$annot_dir
      ~/apps/CAT/CAT_pack/CAT contigs -c ../$ASM \
      -d $CAT_db -t $CATtax_db --force \
      --index_chunks 1 -n $CAT_threads 
      ~/apps/CAT/CAT_pack/CAT add_names -i out.CAT.contig2classification.txt \
      -t $CATtax_db \
      -o CAT.class.out &&\
      mv CAT.class.out \
      $FILE
   fi
done

}



Target_Vir_Contig_Ident () {
echo "
########################################
### Target Vir contig identification ###
############## With BLAST ##############
########################################
"
VIR_REF=$wd/$vir_dir/$vir_target_seq
#make a blast DB for searching with the loop
cd $wd/$vir_dir/
makeblastdb -in $VIR_REF -dbtype nucl
cd $wd
for I in subsmpl-*;
do
   cd $wd/$I
   #Identify target viral sequences in assembly
   #blast assembled contigs against reference viral genomes
   #####Need to add a bit to check that FILE2 is not empty
   func="Identification of assembled target viral contigs in $I"
   ASM=$asm_dir/scaffolds.fasta
   VIR_REF=$VIR_REF
   FILE="${wd}/${I}/${annot_dir}/${vir_target_scaffolds}"
   if [ -s "$FILE" ] ; then
      echo "$FILE already generated (and not empty)
      Skipping identification of assembled target viral contigs
      "
   else
      echo "
      Start ${func}
      $(date)

      Running BLASTN with $VIR_REF as the database and $ASM as query
      "
      blastn -db $VIR_REF -query $ASM \
      -outfmt '6 qaccver saccver pident length qlen mismatch gapopen qstart qend sstart send evalue' \
      > $annot_dir/target_vir.blst6

      #filter blast hits
      cd $wd/$I/$annot_dir
      pident=99
      perc_cov=0.60
      cat target_vir.blst6 | \
      sort -grk 3,3 | awk '{if ($3 >= 98)  print $0}' \
      > target_vir.blst6.filter

      #create file with reference \t contig identity
      cat target_vir.blst6.filter | \
      awk '{print $2,"\t", $1 }' | sort -k 1 \
      > target_vir.blst6.filter.names

      #create file with just contig names for searching
      cat target_vir.blst6.filter.names | awk '{print $2}' | sort\
      > target_vir.blst6.filter.contig.names
      num_target_vir=$(wc -l target_vir.blst6.filter.contig.names)

      #make a fasta file of viral contigs
      ~/apps/bbmap/filterbyname.sh include=T \
      names=target_vir.blst6.filter.contig.names \
      in=$wd/$I/$asm_dir/scaffolds.fasta \
      out=$FILE
   fi
done

}



Ident_Target_Vir_Reads () {
echo "
#########################################################
############# Align reads to assembly then ##############
########## pull out Target viral alignments ############
#########################################################
"
#dowsnt use the unpaired reads
cd $wd
for I in subsmpl-*;
do
   I1=${I#*-}
   I2=${I1%/}
   cd $wd/$I
   #align reads to scaffolds/contigs and compare to Gold standard alignment provided by ART
   #   only compares number
   #  will add a false pos/neg analysis at some point
   #  from that we could maybe check out the character of unassembled viral reads
   func="Aligning Simulated reads to assembly for $I"
   FILE0="${wd}/${I}/${annot_dir}/reads_to_scaf.bam"
   FILE="${wd}/${I}/$reads_dir/${I2}_viral_reads.names"
   #in variables
   REF=$asm_dir/scaffolds.fasta
   READS1=$wd/$I/$reads_dir/bknd-*-vir-*-1.fq
   READS2=$wd/$I/$reads_dir/bknd-*-vir-*-2.fq
   #READS0=$wd/$I/$reads_dir/Up-bknd-*-vir-*-0.fq
   if [ -f "$FILE0" ] ; then
      echo "$FILE0 already generated (and not empty)
      Skipping ${func}
      "
   else
      echo "
      Start ${func}
      $(date)"
      ~/apps/bwa-0.7.17/bwa index $asm_dir/scaffolds.fasta
      ~/apps/bwa-0.7.17/bwa mem -v 1 -t $threads $asm_dir/scaffolds.fasta \
      ${READS1} \
      ${READS2} |\
      samtools sort -n > ${FILE0}.tmp && \
      mv ${FILE0}.tmp $FILE0
   fi


   if [ -f "$FILE" ] ; then
      echo "$FILE and already generated
      Skipping
      "
   else 
      # Filter Art generated reads pull out target viral reads
      # Probably should figure out how to filter out the reads removed during trimming
      cat $reads_dir/${I2}-1.fq  $reads_dir/${I2}-2.fq  |\
      grep -f ../$vir_dir/vir_accession_list.clean |\
      awk -F/ '{print $1}' |\
      sed "s/@//g" |\
      sort | uniq \
      > $FILE #FILE="${wd}/${I}/$reads_dir/${I2}_viral_reads.names"
   fi
done
}



Compare_Vir_Ident_to_Gold_Standard () {
echo "
############################################
####### search annotation files for ########
#### viral contigs identified by blast  ####
############################################
"
#Create final output info
#  Could make it restartable, but this is a very fast module
#  and the whole reason for running the pipe
#pretty lame at the moment, just finds num of target viral seqs correctly identified in DVF and CAT annotations
# With misidentified contigs (target virus called something else) for CAT
out_file=$wd/$vir_ident_summary
echo -e 'file\t#_vir_contigs\tCorrect_ident\tMissclassified\n' >\
$out_file
cd $wd/

for J in $bknd_cov_fractions
do
   for K in $vir_cov_fractions
   do
      I=subsmpl-bknd-$J-vir-$K
      echo "working on $I"
      I1=${I#*-}
      I2=${I1%/}
      cd $wd/$I/
      if [ -s "$PLACE_HOLDER" ]
      then
         echo "#########################"
         echo "$PLACE_HOLDER already preasent"
         echo "#########################"
      else
         cd $wd/$I/$annot_dir
         #TMP
         mv target_vir.blst6.filter.contig.names target_vir.blst6.filter.contig.names.tmp
         cat target_vir.blst6.filter.contig.names.tmp | sort > target_vir.blst6.filter.contig.names
         for L in $DVF_out $VS_out $CAT_out
         do
            cat $L | grep -f target_vir.blst6.filter.contig.names |\
            awk '{print $1}' | sort \
            > ${L%.*}.contig.names
            echo -e\
            "$L\t\
            $J\t\
            $K\t\
            $(cat  target_vir.blst6.filter.contig.names | wc -l)\t\
            $(comm -12 ${L%.*}.contig.names target_vir.blst6.filter.contig.names | wc -l)\t\
            $(comm -13 ${L%.*}.contig.names target_vir.blst6.filter.contig.names | wc -l)" \
            >> $out_file
         done
      fi
   done
done

   echo "##########################"
   echo "DeepVirFinder identified"
   echo "$num_dvf_target of $num_target_vir target viral contigs"
   echo "DVF missed"
   echo "$num_dvf_target_not_found of $num_target_vir target viral contigs"

   ################################################
   ################################################
   #process VirSorter annot
   #could add a score filter
   echo "process VirSorter annot.."
   cat $VS_out | grep -f target_vir.blst6.filter.contig.names |
   awk '{print $1}'\
   > ${VS_out%.*}.contig.names
   num_vs_target=$(cat ${VS_out%.*}.contig.names | wc -l)

   #get viral refs not identified
   cat $VS_out | awk '{print $1}' > ${VS_out%.*}.names
   cat target_vir.blst6.filter.contig.names | grep -vf ${VS_out%.*}.names \
   > ${VS_out%.*}.not_found.names
   num_vs_target_not_found=$(cat ${VS_out%.*}.not_found.names | wc -l)

   echo ''
   echo "##########################"
   echo "VirSorter identified"
   echo "$num_vs_target of $num_target_vir target viral contigs"
   echo "VirSorter missed"
   echo "$num_vs_target_not_found of $num_target_vir target viral contigs"

   ################################################
   ################################################
   #process CAT output
   echo "process CAT annot.."
   cat $CAT_out | grep -f target_vir.blst6.filter.contig.names \
   > ${CAT_out%.*}.target_vir
   cat ${CAT_out}| awk '{print $1}' > ${CAT_out%.*}.all_names
   cat ${CAT_out%.*}.target_vir | awk '{print $1}' > ${CAT_out%.*}.contig.names
   num_CAT_target_virus=$(cat ${CAT_out%.*}.target_vir | grep -i virus | wc -l)
   num_CAT_target_virus_unclassified=$(cat target_vir.blst6.filter.contig.names | grep -vf ${CAT_out%.*}.all_names)
   num_CAT_target_virus_misclassified=$(cat ${CAT_out%.*}.target_vir | grep -iv virus | wc -l)

   num_CAT_target=$(cat ${CAT_out%.*}.target_vir | wc -l)
   echo "#########################"
   echo "CAT identified"
   echo "$num_CAT_target_virus of $num_target_vir target viral contigs"
   echo "CAT mis-identified"
   echo "$num_CAT_target_virus_misclassified of $num_target_vir target viral contigs"
   echo #########################

   #write all stats to summary file $out_file
   # in format 'subsample\t#_vir_contigs\tDeepVirFinder_ident\tDeepVirFinder_missed\tVirSorter_ident\tVirSorter_missed\tCAT_ident\tCAT_unclasifies\tCAT_misclassified\t'

   echo -e \
      "$I\t\
      $num_target_vir\t\
      $num_dvf_target\t\
      $num_dvf_target_not_found\t\
      $num_vs_target\t\
      $num_vs_target_not_found\t\
      $num_CAT_target_virus\t\
      $num_CAT_target_virus_unclassified\t\
      $num_CAT_target_virus_misclassified\t" \
      >> $out_file
}



Pull_Vir_Identified_Contigs () {
echo "
###################################################################
########  Extract viral contigs for each Virus Identifier  ########
###################################################################
"
out_file=$wd/contigs.complete
cd $wd/
for I in subsmpl-*;
do
   out_file=$wd/$contigs_dir/contigs.$I.complete
   func="Extract contigs and place into $wd/$contigs_dir/TOOL.$I.contigs.fasta"
   if [ -f $out_file ]
   then
      echo "Already completed
      $func"
   else
      echo ""
      echo "######################"
      echo "$func"
      echo "######################
      "
      cd $wd/$I/$annot_dir
      #make fasta of target viral contigs
      #have to change names manually...
      cp ../$asm_dir/scaffolds.fasta $wd/$contigs_dir/ALL.$I.contig.fasta && \
      filterbyname.sh in=../$asm_dir/scaffolds.fasta names=target_vir.blst6.filter.contig.names include=true overwrite=true \
      out=$wd/$contigs_dir/BLAST.$I.contig.tmp.fa && mv $wd/$contigs_dir/BLAST.$I.contig.tmp.fa $wd/$contigs_dir/BLAST.$I.contig.fa && \
      filterbyname.sh in=../$asm_dir/scaffolds.fasta names=${DVF_out%.*}.contig.names include=true overwrite=true \
      out=$wd/$contigs_dir/DVF.$I.contig.tmp.fa && mv $wd/$contigs_dir/DVF.$I.contig.tmp.fa $wd/$contigs_dir/DVF.$I.contig.fa && \
      filterbyname.sh in=../$asm_dir/scaffolds.fasta names=${VS_out%.*}.contig.names include=true overwrite=true \
      out=$wd/$contigs_dir/VS.$I.contig.tmp.fa && mv $wd/$contigs_dir/VS.$I.contig.tmp.fa $wd/$contigs_dir/VS.$I.contig.fa && \
      filterbyname.sh in=../$asm_dir/scaffolds.fasta names=${CAT_out%.*}.contig.names include=true overwrite=true \
      out=$wd/$contigs_dir/CAT.$I.contig.tmp.fa && mv $wd/$contigs_dir/CAT.$I.contig.tmp.fa $wd/$contigs_dir/CAT.$I.contig.fa && \
      touch $out_file
   fi
done

}
#############################################



Quantify_Aligned_Reads () {
echo "
##########################################################################
#######  Quantify Reads aligning to contigs from vir Identifiyers  #######
##########################################################################
"

#TODO:this all ignores read pairs...fix it
echo "
Filtering alignments for contigs from different viral identifiers
and make a file of synthetic read names
"

for J in $bknd_cov_fractions 
do
   for K in $vir_cov_fractions
   do
      I=subsmpl-bknd-$J-vir-$K
      echo " working on $I"
      cd $wd/$I/$annot_dir/
      ALIGNMENT="${wd}/${I}/${annot_dir}/reads_to_scaf.bam"
      for L in *.contig.names
      do
         samtools view $ALIGNMENT | grep -f $L |\
         awk '{print $1}' | sort | uniq \
         > ${L%.*.*}.vir_reads
      done
   done
done


echo ""
echo "############################"
echo "Tally up total vir reads, vir reads aligning to vir contigs, vir reads not aligning to vir contigs"
echo "and print to $out_file"
cd $wd
out_file=${wd}/$read_align_stats
touch $out_file
echo -e 'file\tbackground\tvirus\tsim_virus_reads\tcontig_aligning_reads\tunmapped_reads' >\
$out_file
for J in $bknd_cov_fractions
do
   for K in $vir_cov_fractions
   do
      I=subsmpl-bknd-$J-vir-$K
      echo "working on $I"
      I1=${I#*-}
      I2=${I1%/}
      cd $wd/$I/
      if [ -s "$PLACE_HOLDER" ]
      then
         echo "#########################"
         echo "$PLACE_HOLDER already preasent"
         echo "#########################"
      else
         for L in $annot_dir/*.read_names
         do
            echo -e\
            "$L\t\
            $J\t\
            $K\t\
            $(cat  $reads_dir/${I2}_viral_reads.names | wc -l)\t\
            $(comm -12 $L $reads_dir/${I2}_viral_reads.names | wc -l)\t\
            $(comm -13 $L $reads_dir/${I2}_viral_reads.names | wc -l)" \
            >> $out_file 
         done
      fi
   done
done

}
############################################



Calculate_Vir_Contig_Stats () {
echo "
###################################################################
#######  Calculate contigs stats from all vir Identifiyers  #######
###################################################################
"
cd ${wd}/$contigs_dir/
statswrapper.sh ../$vir_dir/$vir_target_seq *.contig.fa minscaf=500 format=4 > vir_contigs_stats.tab 

echo "
############################################
########### Do complete analysis ###########
############## with metaQuast ##############
############################################
"
cd $wd
for I in subsmpl* 
do 
   cd $wd/$I/$asm_dir 
   ln -s scaffolds.fasta $I.fasta 
   cd $wd
done

cd $wd
FILE0=$wd/metaquast_out.single_vir.$run_name
func0="Run mataquast on all viral sequences separately to create $FILE0" 
FILE1=$wd/metaquast_out.all_vir.$run_name
func1="Run mataquast on all viral sequences together to create $FILE1" 


if [ -d "$FILE0" ] ; then
   echo "$FILE0 already generated (and not empty)
   Skipping ${func0}
   "
else
   metaquast.py ./$contigs_dir/*contig.fa -o ${FILE0}.tmp/ -r ./$vir_dir/$single_vir/ -t $threads && \
   mv ${FILE0}.tmp/ $FILE0
fi

# if [ -d "$FILE1" ] ; then
#    echo "$FILE1 already generated (and not empty)
#    Skipping ${func1}
#    "
# else
#    metaquast.py ./$contigs_dir/*contig.fa -o ${FILE1}.tmp/ -r ./$vir_dir/$vir_target_seq -t $threads && \
#    mv ${FILE1}.tmp/ $FILE1
# fi 

}
######################################################

Wrap_Up () {
echo "
######################################################################
################### The full pipeline is finished ####################
############## (brace yourself, it probably didnt work) ##############
######################################################################
"
#ToDO
# add intermediate file cleanup
}



FULL_PIPE
