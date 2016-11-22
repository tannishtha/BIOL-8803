#!/bin/bash
#Here getopts gets the user input
while getopts "n:m:v" opt
do
	case $opt in
   		n)
   		numFiles=$OPTARG
     		;;
   		m)
     		numSeq=$OPTARG
     		;;
   		v)
   	 	verbose=1
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

#The user inputs are used to generate the files 
if [ $numFiles -a $numSeq ]
then
	#Here i contains the number of multi FASTA sequence files given by the user
	for (( i = 1; i <= $numFiles; i++ )) 
     	do
       		a="seq$i.fasta"
       		if [ -f $a ] #here the if block checks if the file is present or not. If it is present, then it is removed.
       		then
	      		rm $a 
	      		if [ $verbose ]
	      		then
	      	  		echo "Deleting file $a as it is already present"
	      		fi
       		fi
       		if [ $verbose ] #New file creation starts from this point
	     	then
	      		echo "Creating new file $a"
	   	fi
	   	#j contains the number of FASTA sequence to be generated in each file
       		for (( j = 1; j <= $numSeq; j++ ))
       		do  
     	 		echo ">seq$i""_$j" >> $a
     	 		if [ $verbose ] 
	     		then
	      	  		echo "Creating fasta sequence >seq$i""_$j"
	     		fi 
     	 		dna=$(cat /dev/urandom | tr -dc 'ACGT' | fold -w 50 | head)  #This generates the random DNA sequence
     	 		echo "$dna" >> $a
       		done
     	done
fi
