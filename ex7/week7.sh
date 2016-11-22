#!/bin/bash

#initialising the variables
presentFile1=0
readFile1=""
presentFile2=0
readFile2=""
presentReference=0
referenceFile=""
outputFilePresent=0
outputFile=""

#Getting the options from the command line
while getopts "a:b:g:o:ezvih" opt
do
	case $opt in
	a)
	   readFile1=$OPTARG
	   presentFile1=1
	   ;;
	b)
	   readFile2=$OPTARG
	   presentFile2=1
	   ;;
	g)
	   referenceFile=$OPTARG
	   presentReference=1
	   ;;
	o)
	   outputFile=$OPTARG
	   outputFilePresent=1
	   ;;
	e)
	   dorealign=1
	   ;;
	z)
	   outzipped=1
	   ;;
	v)
	   verbose=1
	   ;;
	i)
	   indexOutFile=1
	   ;;
	h)
	   echo -e "Usage: ./week7.sh [options]"
	   echo "Options:"
	   echo -e "-a\tInput read file1(required) "
	   echo -e "-b\tInput read file2(required) "
	   echo -e "-g\tReference genome file(required) "
	   echo -e "-o\tOutput VCF file(required) "
	   echo -e "-e\tDo reads re-alignment"
	   echo -e "-z\tGunzip output VCF file"
	   echo -e "-v\tVerbose Mode"
	   echo -e "-i\tIndex output BAM file"
	   echo -e "-h\tPrint usage information" 
	   exit 1
	   ;;
	:) 
   	  echo "Option $OPTARG requires argument" &>2
      exit 1
      ;;
   \?) 
      echo "Invalid option" &>2
      exit 1
      ;;
	esac
done

#to check if all required files are present or not
if [ "$presentFile1" == 1 ] && [ "$presentFile2" == 1 ] && [ "$presentReference" == 1 ] && [ "$outputFilePresent" == 1 ]
then
	if [ -e $readFile1 ] && [ -e $readFile2 ] && [ -e $referenceFile ] 
	then 
		if [ -e $outputFile ]
		then 
			echo "Output file already exists. Do you want to overwrite or exit the program ?(Type overwrite/exit) : "
			read str
			if [ "$str" == "exit" ]
			then
				exit
			fi
		fi
	else
		if [ ! -e $readFile1 ]
		then
			echo "$readFile1 does not exist"
	     	exit
	    elif [ ! -e $readFile2 ]
		then
			echo "$readFile2 does not exist"
			exit
		else
			echo "$referenceFile does not exist"
			exit
		fi
	fi
else
	if [ "$presentFile1" == 0 ]
	then
		echo "Input read file1 not provided. Please use option -a to input it."
		exit
	fi
	if [ "$presentFile2" == 0 ]
	then
		echo "Input read file2 not provided. Please use option -b to input it."
		exit
	fi
	if [ "$presentReference" == 0 ]
	then 
		echo "Reference genome file not provided. Please use option -g to input it."
		exit
	fi
	if [ "$outputFilePresent" == 0 ]
	then 
		echo "Output file name not provided. Please use option -o to input it."
		exit
	fi
	
fi

#Indexing the reference file so as to prepare it for mapping
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nIndexing the reference file.."
fi
bwa index $referenceFile 

#map the reads to the reference
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nMapping reads to reference file to make a SAM file.."
fi
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $referenceFile $readFile1 $readFile2 > lane.sam

#clean read pairing information and flags as BWA can sometimes leave usual flag information on SAM records
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nCleaning read pairing information and flags and storing it as BAM file.."
fi
samtools fixmate -O bam lane.sam lane_fixmate.bam

#sort them from name order to coordinate order
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nSorting in coordinate order.."
fi
samtools sort -O bam -o lane_sorted.bam -T /tmp/lane_temp lane_fixmate.bam

#Relignment of raw gapped alignment
if [ "$dorealign" == 1 ]
then 
	
	# Creating fasta index file for GTAK analysis
	if [ "$verbose" == 1 ]
	then 
		echo -e "\n\nCreating fasta index file.."
	fi
	samtools faidx $referenceFile 
	
	#Creating fasta sequence dictionary file
	if [ "$verbose" == 1 ]
	then 
		echo -e "\n\nCreating fasta sequence dictionary file.."
	fi
	
	temp=$(echo $referenceFile | cut -d"." -f2)
	outFil=".""$temp"".dict"
	java -jar ./dependencies/picard.jar CreateSequenceDictionary REFERENCE=$referenceFile OUTPUT=$outFil
	
	#indexing the bam file
	if [ "$verbose" == 1 ]
	then 
		echo -e "\n\nIndexing the BAM file.."
	fi
	samtools index lane_sorted.bam
	
	if [ "$verbose" == 1 ]
	then 
		echo -e "\n\nPerforming read realignment.."
	fi
	java -Xmx2g -jar ./dependencies/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $referenceFile -I lane_sorted.bam -o lane.intervals --known ./data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	java -Xmx4g -jar ./dependencies/GenomeAnalysisTK.jar -T IndelRealigner -R $referenceFile -I lane_sorted.bam -targetIntervals lane.intervals -known ./data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o lane_realigned.bam
fi

#Index the output BAM file
if [ "$indexOutFile" == 1 ]
then
	if [ "$verbose" == 1 ]
	then 
		echo -e "\n\nIndexing output BAM file.."
	fi
	if [ "$dorealign" == 1 ]
	then
		samtools index lane_realigned.bam
	else
		samtools index lane_sorted.bam
	fi
fi

#Variant calling
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nVariant calling.."
fi
if [ "$indexOutFile" == 1 ] && [ "$dorealign" == 1 ]
then 
	samtools mpileup -go study.bcf -f $referenceFile lane_realigned.bam
elif [ "$indexOutFile" == 1 ]
then
	samtools mpileup -go study.bcf -f $referenceFile lane_sorted.bam
else
	samtools mpileup -go study.bcf -f $referenceFile lane_sorted.bam
fi
if [ "$outzipped" == 1 ]
then
	outputFile="$outputFile"".gz"
	bcftools call -vmO z -o $outputFile study.bcf #z is for compressed VCF
else
	bcftools call -vmO v -o $outputFile study.bcf #v is for uncompressed VCF
fi

#To prepare our VCF for querying we next index it using tabix:
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nIndexing the VCF file.."
fi
tabix -p vcf $outputFile

#Filter VCF file 
if [ "$verbose" == 1 ]
then 
	echo -e "\n\nFiltering the VCF file.."
fi
bcftools filter -O z -o outputFilter.vcf.gz -s LOWQUAL -i'%QUAL>10' $outputFile


