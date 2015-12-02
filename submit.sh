#!/bin/bash

#Name of Job - will appear in qstat
#$ -N AWS-transfer

#Join your STDOUT and STDERR reports (if processing 50+ files this option is handy but not required)
#$ -j y

#Path of executable (in this case Perl)
#$ -S /usr/bin/perl

#Use current working directory
#$ -cwd

#Export environment variables to qsub
#$ -V

#Number of queues (CPUS) to use, batch size
#$ -t 1

#Number of files in batch size to submit at once
#$ -tc 1

#Number of threads to use per CPU (50 of the 52 nodes have 16 threads available)
#none

#sumbit to local submitting node (for spawning more jobs)
#not needed for this job: -l h=ihg-node-27


#defining SGE_INDEX and input files
#SAMPLE_DIR=/mnt/state_lab/proc/REI1344_SFARI/combinedVCF/combine_var; #define this directory for your files
#SAMPLE_LIST=($SAMPLE_DIR/*.bam);
#INDEX=$((SGE_TASK_ID-1));
#INPUT_FILE=${SAMPLE_LIST[$INDEX]};

#Script to execute
perl /home/smithl/pipeline/aws-transfer/aws_transfer_WGS.pl
