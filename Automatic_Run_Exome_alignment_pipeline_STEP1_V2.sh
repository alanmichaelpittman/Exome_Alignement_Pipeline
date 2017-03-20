#!/bin/bash

#pipeline to align WGS data from HiSeq FastQ files
#Alan Pittman October 2016

#Automatically run Exome_Aligment_pipeline
########################################################################################

#1) This script must be copied into and run (sh) in the run folder

#2) INSERT PROJECT FOLDER HERE: REMEMBER TO INCLUDE DATE !!! TAKE THE RUN DATE EG - 20161012

#3) Also helpful to rename the script with the project being analysed



while getopts p: option
do
        case "${option}"
        in
                p) PROJECT=${OPTARG};;


  esac

done


########################################################################################

#get input and output directories

MasterDIR=`pwd`

echo "MASTER DIRECTORY:"

echo " "
echo $MasterDIR
echo " "

inputDIR="$MasterDIR/Unaligned/$PROJECT"
outputDIR="$MasterDIR/a_aligned_X/$PROJECT"

mkdir $MasterDIR/a_aligned_X

echo "PROJECT DIRECTORY:"
echo " "
echo $inputDIR
echo " "
echo $outputDIR
echo " "

#get sample ID's

MySampleIDs=`ls $MasterDIR/Unaligned/$PROJECT`

echo "SAMPLE IDs:"
echo $MySampleIDs
echo " "

echo "PIPELINE!"

sleep 5

#and now finally a for loop with all the above variables to invoke the master script:

for nID in $MySampleIDs; do
	sh /array/HiSeq3000/Pipeline_scripts/Exome_alignment_pipeline_STEP1.sh -i $inputDIR -o $outputDIR -l $nID
done

exit

########################################################################################
























novoalign="/data/kronos/NGS_Software/novocraft_v3/novoalign"
novoalignArguments="-c 8 --rOQ --hdrhd 3 -H -k -o Soft -t 320 -F ILM1.8"
indexedgenome="/data/kronos/NGS_Reference/novoalign/human_g1k_v37.fasta.k15.s2.novoindex"
samtools="/data/kronos/NGS_Software/samtools-0.1.18/samtools"
picard="/data/kronos/NGS_Software/picard-tools-1.75"
gatk="/data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar"
genomeFASTA="/data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta"
genomeFASTAI="/data/kronos/NGS_Reference/fasta/human_g1k_v37.fastai"
JAVA="/data/kronos/General_Software/jre1.7.0_67/bin/java"
knownINDELS="/data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf"
samtoolsMEM="5000000000" 
picardArguments="TMP_DIR=/data/kronos/temp ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE" 

########################################################################################

echo `date`

echo "
	====================================
	Whole Genome Sequencing Alignment 
	====================================
	====================================
	Alan Pittman - November 2015 
	====================================
	"	

sleep 1

echo "
	===================================
	STEP ONE - checking for fastq files
	===================================
	"
sleep 2
	
########################################################################################

while getopts i:o:l: option
do
        case "${option}"
        in
                i) iFolder=${OPTARG};;
                o) oFolder=${OPTARG};;
				l) myIDs=${OPTARG};;

  esac

done

if [ ! -e $oFolder ]; then mkdir $oFolder; echo "making output directory"; fi

########################################################################################

n=0
while (( $n <= 7 ))
	do
	n=$(( n+1 ))

if [ -e ${iFolder}/${myIDs}/*L00${n}_R1*gz ]; then
	
	echo "Files exist"
	echo "The are FASTQ files from lane ${n}"

	for nID in $myIDs; do

		folder=${iFolder}/${nID}
		
		for file1 in `find $folder -name *L00${n}_R1*gz`; do
		((nfiles=nfiles+2))
		file2=`echo $file1 | sed -e 's/R1/R2/g'`
		inputFiles=`echo $file1 $file2`
		
		echo $inputFiles
		
		output=${oFolder}/${nID}/${nID}
		
		if [ ! -e ${oFolder}/${nID} ]; then mkdir ${oFolder}/${nID}; fi 


sleep 3	

		
echo "
	===================================
	STEP TWO - writing the scripts
	===================================
	"
sleep 1

echo "
	done !
	"
sleep 1				
		echo "#!/bin/bash"  > $oFolder/J_script_${nID}_${n}.sh
		
		echo "#These are the pipeline steps for sample ${nID} of lane ${n} of the HiSeq"  >> $oFolder/J_script_${nID}_${n}.sh
			
##now lets write to script each step in the pipeline:		
		
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$novoalign -o SAM $'@RG\tID:${nID}\tSM:${nID}\tLB:${nID}\tPL:ILLUMINA' -f $inputFiles -d $indexedgenome $novoalignArguments > ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools view -bS -t $genomeFASTAI -o ${output}_${n}.bam ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "rm ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools sort -m $samtoolsMEM ${output}_${n}.bam ${output}_${n}_sorted" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools index ${output}_${n}_sorted.bam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "java -Xmx10g -jar $picard/MarkDuplicates.jar $picardArguments INPUT=${output}_${n}_sorted.bam OUTPUT=${output}_${n}_sorted_unique.bam METRICS_FILE=${output}_${n}_picard_metrics.out" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools index ${output}_${n}_sorted_unique.bam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo $JAVA -jar $gatk -T RealignerTargetCreator	-nt 8 -R $genomeFASTA -o ${output}_${n}_sorted_unique.bam.list -I ${output}_${n}_sorted_unique.bam --known $knownINDELS >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo $JAVA -jar $gatk -T IndelRealigner -R $genomeFASTA -targetIntervals ${output}_${n}_sorted_unique.bam.list -I ${output}_${n}_sorted_unique.bam -o ${output}_${n}_sorted_unique_realigned.bam --knownAlleles $knownINDELS >> $oFolder/J_script_${nID}_${n}.sh
	
##removing intermediate files:	
	
		echo "
	
			###Remove Intermediate BAM files leaving only final processed_GATK_bam file:

			if [ -e ${output}_${n}_sorted_unique_realigned.bam ] 
				then
					rm ${output}_${n}.bam
					rm ${output}_${n}.sam
					rm ${output}_${n}_sorted.bam
					rm ${output}_${n}_sorted.bam.bai
					rm ${output}_${n}_sorted_unique.bam
					rm ${output}_${n}_sorted_unique.bam.bai
					rm ${output}_${n}_sorted_unique.bam.list
					
					
					echo 'old BAM and intermediate BAM.bai all removed'
				else
					echo '! something has gone wrong '
				fi

		

			" >> $oFolder/J_script_${nID}_${n}.sh
			
			
			echo "exit" >> $oFolder/J_script_${nID}_${n}.sh
	
		done
		
	done
		
	else echo "There are NO FASTQ files from lane ${n}"
	sleep 1
fi	

done


##now last (but not least) we will automatically submit all jobs to the SGE of KRONOS

	echo "
	=======================================
	STEP THREE - Submitting jobs to Kronons
	=======================================
	"
sleep 1
	
n=0
while (( $n <= 7 ))
	do
	n=$(( n+1 ))

		for nID in $myIDs; do

			if [ -e ${iFolder}/${myIDs}/*L00${n}_R1*gz ]; then
				sleep 2
				qsub -pe make 18 -cwd $oFolder/J_script_${nID}_${n}.sh
				
			else 
			sleep 2
			echo "There is NO job from lane ${n}"

			fi
		done	

done

exit



#misc.commands:


-o SAM $'@RG\tID:hmorris-IONNG-79658\tSM:hmorris-IONNG-79658\tLB:hmorris-IONNG-79658\tPL:ILLUMINA' 


/data/kronos/General_Software/jre1.7.0_67/bin/java -jar /data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar -T RealignerTargetCreator -L /data/kronos/NGS_Reference/GATK_refFiles/INTERVALS.bed  -nt 5 -R /data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta -o /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam.list -I /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam --known /data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf 	
/data/kronos/General_Software/jre1.7.0_67/bin/java -jar /data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar -T IndelRealigner -targetIntervals /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam.list -I /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam  -R /data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta -o /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique_realigned.bam --knownAlleles /data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf 
java -Xmx10g -jar /data/kronos/NGS_Software/picard-tools-1.75/MarkDuplicates.jar 
TMP_DIR=/data/kronos/temp ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE 
INPUT=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted.bam 
OUTPUT=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam 
METRICS_FILE=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_picard_metrics.out
novoalign arguments 
-c 8 --rOQ --hdrhd 3 -H -k -o Soft -t 320 -F ILM1.8 -d /data/kronos/NGS_Reference/novoalign/human_g1k_v37.fasta.k15.s2.novoindex
/data/kronos/NGS_Software/novocraft_v3/novoalign -c 8 -o SAM $'@RG\tID:2PhCl\tSM:2PhCl\tLB:2PhCl\tPL:ILLUMINA' --rOQ --hdrhd 3 -H -k -a CTGTCTCTTATACACATCT -o Soft -t 320 -F ILM1.8 -f /data/kronos/apittman/Exome_Data_Analysis/Christos/Unaligned_Exp2/2PhCl/2PhCl_S2_L008_R1_001.fastq.gz /data/kronos/apittman/Exome_Data_Analysis/Christos/Unaligned_Exp2/2PhCl/2PhCl_S2_L008_R2_001.fastq.gz  -d /data/kronos/NGS_Reference/novoalign/human_g1k_v37.fasta.k15.s2.novoindex > /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl.sam
/data/kronos/NGS_Software/samtools-0.1.18/samtools view -bS -t /data/kronos/NGS_Reference/fasta/human_g1k_v37.fastai -o /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl.bam /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl.sam  ## make BAM file
/data/kronos/NGS_Software/samtools-0.1.18/samtools sort -m 5000000000 /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl.bam /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted  ## sort
/data/kronos/NGS_Software/samtools-0.1.18/samtools index  /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted.bam   ##build index
## Now remove PCR duplicates using PICARD
java -Xmx10g -jar /data/kronos/NGS_Software/picard-tools-1.75/MarkDuplicates.jar TMP_DIR=/data/kronos/temp ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted.bam OUTPUT=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam METRICS_FILE=/data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_picard_metrics.out
/data/kronos/NGS_Software/samtools-0.1.18/samtools index  /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam 
#Alan Pittman EXTRA STEPS TO PERFORM IndelREALIGNMENT AND VARIANT DETECTION IN GATK
/data/kronos/General_Software/jre1.7.0_67/bin/java -jar /data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar -T RealignerTargetCreator -L /data/kronos/NGS_Reference/GATK_refFiles/INTERVALS.bed  -nt 5 -R /data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta -o /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam.list -I /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam --known /data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf 	
/data/kronos/General_Software/jre1.7.0_67/bin/java -jar /data/kronos/NGS_Software/GATK_v3_3/GenomeAnalysisTK.jar -T IndelRealigner -targetIntervals /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam.list -I /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique.bam  -R /data/kronos/NGS_Reference/fasta/human_g1k_v37.fasta -o /data/kronos/apittman/Exome_Data_Analysis/Christos/a_aligned_v4_Exp2/2PhCl/2PhCl_sorted_unique_realigned.bam --knownAlleles /data/kronos/NGS_Reference/GATK_refFiles/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf 

