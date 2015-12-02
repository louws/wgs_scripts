#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Find;
use Data::Dumper; # for debugging

sub help();
sub getPipelineChoice();
sub getHeaderFromFq($$);
sub mergeBamsAndRealign($$);
sub combineGVCFs($$);
sub getSampleIDFromBam($);
sub checkForExistingDirs(@);
sub getLinesForReport($);

## program and reference data paths (to be cleaned up later)
my $picardLoc = "/mnt/state_lab/progs/picard-tools-1.113";
my $refGenome = '/mnt/state_lab/scripts/resources/Pipeline/human_g1k_v37.fasta';
my $bwaBin = '/mnt/state_lab/progs/bwa-0.7.6a/bwa';
my $javaLoc= '/mnt/state_lab/progs/jdk1.7.0_55/bin';
my $goldenIndels = '/mnt/state_lab/reference/recal_training/Mills_and_1000G_gold_standard.indels.b37.vcf';
my $dbSNP138 = '/mnt/state_lab/reference/recal_training/dbsnp_138.b37.vcf';
my $gatk = '/mnt/state_lab/progs/gatk/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar';
my $samtoolsLoc = '/mnt/state_lab/progs/samtools-0.1.19/samtools';
my $recalTraining1 = '/mnt/state_lab/reference/recal_training/hapmap_3.3.b37.vcf';
my $recalTraining2 = '/mnt/state_lab/reference/recal_training/1000G_omni2.5.b37.vcf';
my $recalTraining3 = '/mnt/state_lab/reference/recal_training/1000G_phase1.snps.high_confidence.b37.vcf';
my $recalTraining4 = $dbSNP138;
my $plinkSeq = '/mnt/state_lab/progs/plinkseq/plinkseq-0.10/pseq';
my $pseqGenome = '/mnt/state_lab/proc/willseyJ/plinkseq/hg19';
my $nonPARIntervals = '/mnt/state_lab/scripts/resources/GRCh37_nonPARIntervals.intervals';

# dividing into 13 approximately equivalent regions and preparing appropriate file endings
my @chrGroups = ("'-L 1'", "'-L 2'", "'-L 3 -L 21'", "'-L 4 -L 22'", "'-L 5 -L 19'", "'-L 6 -L Y'", "'-L 7 -L 20 -L MT'", 
	"'-L X -L 18'", "'-L 8 -L 17'", "'-L 9 -L 16'", "'-L 10 -L 15'", "'-L 11 -L 14'", "'-L 12 -L 13'");
my $numChrGroups = scalar @chrGroups;
my @vcfEndings; # will look like this: chr1.vcf, ..., chr8,17.vcf, etc.
for (@chrGroups) {
	my $vcfName = "'chr";
	(my $strip = $_) =~ s/'|-L\s//g;
	my @chrNums = split /\s+/, $strip;
	$vcfName .= (join ",", @chrNums) . ".vcf'";
	push @vcfEndings, $vcfName;
}
my $numJointJobs = scalar @chrGroups;

# adding required items to environment in case they're not there
$ENV{'PATH'} = "$javaLoc:$ENV{'PATH'}"; # gets used internally by GATK (maybe)
$ENV{'LD_LIBRARY_PATH'} = "/opt/sge625/sge/lib/lx24-amd64/:/mnt/state_lab/progs/hdf5/hdf5-1.8.13/hdf5/lib/:/mnt/state_lab/progs/boostlib/lib"; # which is actually used, if any?

# getting options
my ($outDir, $inputDir, $sampleFile, $pedFile, $phenFile, $getHelp);
my $intervals = '';
my $padding = 100; # default value
my ($fqForParsing, $mergeSameSampleBams, $combineGVCFs, $bamToCheckID);
my $maxJobs = 91; # default value for qsub "tc" option

my %options = (
	# always required
	'out:s' => \$outDir,
	'in:s' => \$inputDir,

	# required for specific pipeline steps
	'sample-list|sampleList|list:s' => \$sampleFile,
	'pedigree:s' => \$pedFile,
	'phenotype|phe:s' => \$phenFile, # for PlinkSeq

	# optional
	'intervals:s' => \$intervals,
	#'padding:s' => \$padding, (currently letting it always be 100)
	'max-jobs|maxjobs:s' => \$maxJobs,
	'h|help|u' => \$getHelp,

	# for internal use
	'getHeaderFromFq:s' => \$fqForParsing,
	'mergeSameSampleBams' => \$mergeSameSampleBams,
	'combineGVCFs' => \$combineGVCFs,
	'getSampleIDFromBam:s' => \$bamToCheckID,
);
GetOptions(%options);

$0 = abs_path $0; # this is actually necessary in some cases!

# these subroutines run when the program calls itself midway through the pipeline
# for the most part, the arguments ($inputDir, $outDir, etc.) should already be validated
if ($fqForParsing) {
	getHeaderFromFq($fqForParsing, $sampleFile);
	exit;
} elsif ($mergeSameSampleBams) {
	mergeBamsAndRealign($inputDir, $outDir);
	exit;
} elsif ($combineGVCFs) {
	combineGVCFs($inputDir, $outDir);
	exit;
} elsif ($bamToCheckID) {
	my $sampleID = getSampleIDFromBam($bamToCheckID);
	print $sampleID;
	exit;
}

# giving help if requested or if input directory is missing
if ($getHelp || !$inputDir) {
	help();
	exit;
}

# if no output directory specified, asking for permission to use current directory
unless ($outDir) {
	print "You haven't specified an output directory with the \"-out\" option.\nIs it okay to place output files inside your current directory (yes/no)? ";
	while () {
		chomp (my $answer = <STDIN>);
		if ($answer =~ /^n(o)?$/i) {
			print "Exiting. No jobs were submitted.\n";
			exit;
		}
		elsif ($answer =~ /^y(es)?$/i)  {
			$outDir = getcwd();
			last;
		} else {
			print "Invalid response. Answer \"yes\" or \"no:\"";
		}
	}
}

# setting up the intervals (-L) option for gatk
my @intervalFiles = split /,/, $intervals;
for (@intervalFiles) {
	if (-e $_) {
		$_ = '-L ' . abs_path $_;
	} else {
		print "Error: Couldn't find interval file $_.\n";
		exit;
	}
}
$intervals = join " ", @intervalFiles;
die "Invalid padding value (must be a number)" unless ($padding =~ /^\d+$/);
$padding = "-ip $padding";


# making sure program is running from a submit node (required for all the job submissions!)
chomp(my $node = `hostname`);
unless ($node =~/^ihg-node-1$|^ihg-node-27$|^ihg-node-50$|^ihg-headnode-1$/) {
	print "This program must run from a node with the power to qsub (nodes 1, 27, or 50, or the headnode).\n";
	print "Log into one of these nodes and run the program again.\n";
	print "Hint: From the head node, the command 'qlogin -l h=ihg-node-50' will log you into node 50.\n";
	exit;
}

# making sure the input/output directory paths are okay
$outDir =~ s/~/$ENV{'HOME'}/; # because the abs_path function doesn't interpret the '~' symbol
$inputDir =~ s/~/$ENV{'HOME'}/;
$inputDir = abs_path $inputDir;
unless (-d $inputDir) {
	print "Invalid input directory; it doesn't seem to exist.\n";
	exit;
}
if (-f $outDir) {
	print "Error: There's an existing file with the same name as your specified output directory.\n";
	exit;
}
mkdir $outDir; # okay if it already exists
$outDir = abs_path($outDir) or die "Couldn't find absolute path of $outDir";
unless (-w $outDir) {
	print "Error: You do not seem to have write permissions for your chosen output directory. Please choose another.\n";
	exit;
}
if ($inputDir eq $outDir) {
	print "Error: The input and output directories cannot be the same directory. You can specifiy any other writable directory,\n";
	print "or supply a path for a new directory to be created.\n";
	exit;
}
my $initialInputDir = $inputDir; # keeping track of this for later

# getting which steps in the pipeline to run from the user
my %stepsToRun = getPipelineChoice();
my $numFiles; # used throughout pipeline to determine number of tasks to run in each array job

my @directoriesToMake; # holds list of directories to make (typically, 2-3 needed per step of pipeline)
my @jobScripts; # holds all of the job scripts, which are submitted at the end

# Step 1: Merge split files, gunzip gzipped files, make symbolic links, make all file extensions .fastq
my $interleaved = -9; # keeping track of whether Fastq files are interleaved or paired-end
if ($stepsToRun{'1'}) {
	my $fqDir = $outDir . '/prepped_fastq';
	my $logDir = $fqDir . '/logs';
	push @directoriesToMake, ($fqDir, $logDir);
	my @files = split /\n/, `find $inputDir -maxdepth 1`;

	my @mergedFiles;
	my %filesToProcess; # keys will be file basenames, values will be arrays of the names of the split files associated with each basename
	for (@files) {
		next unless ($_ =~ /\.(fastq|fq)(\.gz)?$/); # skipping all non-Fastq files
		if ($_ !~ /R[12](_\d{3})?\.(fastq|fq)(\.gz)?$/) {
			if ($interleaved == 0) {
				print "Error: There seems to be a mix of interleaved and paired-end .fastq files in the input directory.\n";
				exit;
			} else {
				$interleaved = 1;
			}
		} else {
			if ($interleaved == 1) {
				print "Error: There seems to be a mix of interleaved and paired-end .fastq files in the input directory.\n";
				exit;
			} else {
				$interleaved = 0;
			}
		}
		(my $basename = $_) =~ s/.*\///; # stripping paths
		(my $rootName = $basename) =~ s/(_\d{3})?\.(fastq|fq)(\.gz)?$//; # stripping extensions and split fastq numbering
		
		push @{$filesToProcess{$rootName}}, $basename;
	}
	if (scalar keys %filesToProcess == 0) {
		print "Error: There do not seem to be any fastq files in the input directory.\n";
		exit;
	}

	my @fqPrepCommands;
	for (keys %filesToProcess) {
		# retrieving the list of files to cat
		my @fileArray = @{$filesToProcess{$_}};
		sub getOrder($) {
			my $file = shift @_;
			(my $num = $file) =~ s/.*_(\d{3})\.(fastq|fq)(\.gz)?$/$1/; # getting just the split fastq numbering
			return $num;
		}
		@fileArray = sort {getOrder($a) <=> getOrder($b)} @fileArray; # crucial because the Unix find command doesn't return a sorted list
		my $fileList = join " ", @fileArray;

		# now creating the name of the final merged file and assembling the cat command
		my $rootName = $_;
		my $mergedFile = $fqDir . '/' . $rootName . '.fastq';
		$mergedFile .= '.gz' if ($fileArray[0] =~ /\.gz$/);
		push @mergedFiles, $mergedFile; # just keeping a count of how many total merged files there will be (array ordering doesn't matter)
		
		if (scalar @fileArray == 1) {
			my $file = shift @fileArray;
			push @fqPrepCommands, "'ln -sT $inputDir/$file $mergedFile'";
			#``; # if file doesn't need merging, will just make a symbolic link
			#exit;
		} else {
			push @fqPrepCommands, "'cat $fileList > $mergedFile'";
		}
	}
	$numFiles = scalar @mergedFiles;
	if ($numFiles % 2 != 0 && $interleaved == 0) {
		print "Error: There appears to be an odd number of paired-end fastq files in the input directory.\n";
		exit;
	}
	$numFiles /= 2 unless ($interleaved); # now recording the number of file pairs for paired-end fastq
	my $numCommands = scalar @fqPrepCommands;
	my $prepCommandList = '(' . (join ' ', @fqPrepCommands) . ')';

	my $runCat = << "	EOF";
		\#\$ -N prepFq -j y -o $logDir/prep.\$TASK_ID.txt -t 1-$numCommands -tc $maxJobs -V -S /bin/bash
		cd $inputDir
		INDEX=\$((SGE_TASK_ID - 1))
		COMMANDS=$prepCommandList
		CURR_COMMAND=\${COMMANDS[\$INDEX]}
		THIS_LOG+=`echo "\$CURR_COMMAND" | sed -e "s/.*\\///"`
		THIS_LOG=\${THIS_LOG/%.gz/}
		THIS_LOG=\${THIS_LOG/%fq/txt}
		THIS_LOG=\${THIS_LOG/%fastq/txt}
		mv "$logDir/prep.\$SGE_TASK_ID.txt" "$logDir/prep.\$SGE_TASK_ID.\$THIS_LOG"
		echo "Running this: \$CURR_COMMAND"
		eval \$CURR_COMMAND
		echo "Done."
	EOF
	$runCat =~ s/\t//g;
	push @jobScripts, $runCat;
	$inputDir = $fqDir;
}

if ($stepsToRun{'2'}) {
	# making sure files end in fastq
	unless ($stepsToRun{'1'}) {
		my @fileCheck = split /\n/, `ls $inputDir/*.fastq* 2> /dev/null`;
		$numFiles = scalar @fileCheck;
		if ($numFiles == 0) {
			print "Error: No files containing \".fastq\" found in directory $inputDir.\n";
			print "If you are starting the pipeline with BWA, all input files must contain the extension \".fastq.\"\n";
			exit;
		}
		if ($fileCheck[0] =~ /_R[12]\.fastq/) { # technically should check all the files
			$interleaved = 0;
		} else {
			$interleaved = 1;
		}
		if ($numFiles % 2 != 0 && !$interleaved) {
			print "Error: There appears to be an odd number of paired-end Fastq files in the input directory.\n";
			exit;
		}
		$numFiles /= 2 unless ($interleaved); # recording number of file pairs for paired-end files
	}

	# Checking the SampleID file for errors (it actually gets read in again later)
	my %sampleIDs;
	unless (-f $sampleFile) {
		print "Error: To work with Fastq files, you must supply a file listing SampleIDs.\n";
		print "Run \"perl $0 -help\" for more information.\n";
		exit;
	} else {
		$sampleFile = abs_path $sampleFile;
	}
	open my $IN, '<', $sampleFile or die $!;
	while (my $line = <$IN>) {
		chomp $line;
		next unless ($line);
		my @fields = split /\t/, $line;
		if (scalar @fields != 2) {
			print "Error: Wrong number of columns in line $. of $sampleFile.\n";
			exit;
		}
		my $file = $fields[0];
		my $sampleID = $fields[1];
		if (exists $sampleIDs{$file}) {
			print "Error: File $file appears more than once in the Fastq-SampleID file.\n";
			exit;
		} else {
			$file =~ s/\.fq(\.gz)?$/.fastq/; # by this point in the pipeline, all files should end in .fastq, so fixing sample list entries that end in .fq
			$sampleIDs{$file} = $sampleID if ($file =~ /_R1(_\d{3})?\.fastq(\.gz)?$/ || $interleaved); # for paired end, only taking the R1 files
		}
	}
	unless (scalar keys %sampleIDs > 0) {
		print "Error: There were no valid Fastq filenames in the Fastq-SampleID file.\n";
		print "Run \"perl $0 -help\" for more information.\n";
		exit;
	}
	close $IN;
	$sampleFile = abs_path $sampleFile; # when the file is used again in another shell, will need the correct absolute path

	# Running BWA (Note: avoiding node 28 for now due to a problem with /usr/bin/lesspipe.sh)
	my $samDir = $outDir . '/sam';
	my $logDir = $samDir . '/logs';
	push @directoriesToMake, ($samDir, $logDir);
	(my $genomeBase = $refGenome) =~ s/\.fasta//; # bwa doesn't want the extension
	my $threads = 3; # 3-5 often best, per Louw
	my $maxBwaJobs = int ($maxJobs / $threads); # need smaller job limit to keep the number of slots the same
	my $runBwa = << "	EOF";
		\#\$ -N bwa -j y -o $logDir/bwa.\$TASK_ID.txt -V -S /bin/bash -t 1-$numFiles -tc $maxBwaJobs -l mem_free=8G -l h=!ihg-node-28 -pe parallel $threads
		SAMPLE_LIST_R1=($inputDir/*_R1*fastq*)
		SAMPLE_LIST_R2=($inputDir/*_R2*fastq*)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE_R1=\${SAMPLE_LIST_R1[\$INDEX]}
		INPUT_FILE_R2=\${SAMPLE_LIST_R2[\$INDEX]}
		BASENAME=`basename \$INPUT_FILE_R1`
		OUTPUT_FILE=`echo \$BASENAME | sed -e "s/_R1.*/.sam/"`
		THIS_LOG=\${OUTPUT_FILE/%sam/txt}
		mv "$logDir/bwa.\$SGE_TASK_ID.txt" "$logDir/bwa.\$SGE_TASK_ID.\$THIS_LOG"
		HEADER=`perl $0 -getHeaderFromFq \$INPUT_FILE_R1 -list $sampleFile`
		echo "Processing \$INPUT_FILE_R1 and \$INPUT_FILE_R2"
		echo "Executing on `hostname`"
		if [ "\$HEADER" == 'ERROR' ]
		then
			echo "There was an error creating the read group header line for BWA."
			echo "Most likely, the input file is missing from the Fastq-SampleID list file."
			exit
		fi
		echo "SAM header will be \$HEADER"
		##REPORT
		#BWA:
		COMMAND="time $bwaBin mem -M -t $threads -R \$HEADER $genomeBase \$INPUT_FILE_R1 \$INPUT_FILE_R2 > $samDir/\$OUTPUT_FILE"
		##END_REPORT
		eval \$COMMAND
	EOF

	# different script for interleaved fastq files
	if ($interleaved) {
		$runBwa = << "		EOF";
			\#\$ -N bwa -j y -o $logDir/bwa.\$TASK_ID.txt -V -S /bin/bash -t 1-$numFiles -tc $maxBwaJobs -l mem_free=8G -l h=!ihg-node-28 -pe parallel $threads
			SAMPLE_LIST=($inputDir/*fastq*)
			INDEX=\$((SGE_TASK_ID-1))
			INPUT_FILE=\${SAMPLE_LIST[\$INDEX]}
			BASENAME=`basename \$INPUT_FILE`
			OUTPUT_FILE=`echo \$BASENAME | sed -r -e "s/\\.fastq(\\.gz)?\$/.sam/"`
			THIS_LOG=\${OUTPUT_FILE/%sam/txt}
			mv "$logDir/bwa.\$SGE_TASK_ID.txt" "$logDir/bwa.\$SGE_TASK_ID.\$THIS_LOG"
			HEADER=`perl $0 -getHeaderFromFq \$INPUT_FILE -list $sampleFile`
			echo "Processing \$INPUT_FILE"
			echo "Executing on `hostname`"
			echo "SAM header will be \$HEADER"
			##REPORT
			#BWA:
			COMMAND="time $bwaBin mem -M -p -t $threads -R \$HEADER $genomeBase \$INPUT_FILE > $samDir/\$OUTPUT_FILE"
			##END_REPORT
			eval \$COMMAND
		EOF
	}
	$runBwa =~ s/\t//g;
	push @jobScripts, $runBwa;
	$inputDir = $samDir; # for the next step in the pipeline
}


# Step 3: Converting SAM files to sorted BAM files
if ($stepsToRun{'3'}) {
	unless ($stepsToRun{'2'}) {
		my @fileCheck = split /\n/, `ls $inputDir/*sam 2> /dev/null`;
		$numFiles = scalar @fileCheck;
		if ($numFiles == 0) {
			print "Error: No files ending in \".sam\" found in directory $inputDir.\n";
			exit;
		}
	}

	my $bamDir = $outDir . '/initial_bam';
	my $logDir = $bamDir . '/logs';
	push @directoriesToMake, ($bamDir, $logDir);
	my $makeBams = << "	EOF";
		\#\$ -N make_bam -j y -o $logDir/sort.\$TASK_ID.txt -V -S /bin/bash -l mem_free=8G -t 1-$numFiles -tc $maxJobs
		cd $inputDir # just makes things easier
		SAM_FILES=(*sam)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${SAM_FILES[\$INDEX]}
		OUTPUT_FILE=`echo \$INPUT_FILE | sed -e "s/\\.sam/.bam/"`
		THIS_LOG=\${OUTPUT_FILE/%bam/txt}
		mv "$logDir/sort.\$SGE_TASK_ID.txt" "$logDir/sort.\$SGE_TASK_ID.\$THIS_LOG"
		JOB_NUM=\$((SGE_TASK_ID))
		WORK_DIR="$bamDir/work/\$JOB_NUM"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		echo "Executing SortSam on `hostname`"
		##REPORT
		#Sort SAM and output as BAM:
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $picardLoc/SortSam.jar INPUT=$inputDir/\$INPUT_FILE OUTPUT=$bamDir/\$OUTPUT_FILE SO=coordinate
	EOF
	$makeBams =~ s/\t//g;
	push @jobScripts, $makeBams;
	$inputDir = $bamDir;
}

# Step 4: Mark duplicates and create BAM index files
if ($stepsToRun{'4'}) {
	unless ($stepsToRun{'3'}) {
		my @fileCheck = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
		$numFiles = scalar @fileCheck;
		if ($numFiles == 0) {
			print "Error: No files ending in \".bam\" found in directory $inputDir.\n";
			exit;
		}
	}

	my $indexedBamDir = $outDir . '/indexed_bam';
	my $logDir = $indexedBamDir .'/logs';
	my $metricsDir = $indexedBamDir . '/metrics';
	push @directoriesToMake, ($indexedBamDir, $logDir, $metricsDir);

	my $indexBams = << "	EOF";
		\#\$ -N index_bam -j y -o $logDir/index.\$TASK_ID.txt -V -S /bin/bash -l mem_free=8G -t 1-$numFiles -tc $maxJobs
		cd $inputDir
		SORTED_FILES=(*bam)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${SORTED_FILES[\$INDEX]}
		THIS_LOG=\${INPUT_FILE/%bam/txt}
		mv "$logDir/index.\$SGE_TASK_ID.txt" "$logDir/index.\$SGE_TASK_ID.\$THIS_LOG"
		OUTPUT_FILE=`echo \$INPUT_FILE | sed -e "s/\\.bam/_indexed.bam/"`
		METRICS=`echo \$INPUT_FILE | sed -e "s/\\.bam/_metrics.txt/"`
		JOB_NUM=\$((SGE_TASK_ID))
		WORK_DIR="$indexedBamDir/work/\$JOB_NUM"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		echo "Executing MarkDuplicates on `hostname`"
		##REPORT
		#MarkDuplicates (and make BAM index):
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $picardLoc/MarkDuplicates.jar METRICS_FILE=$metricsDir/\$METRICS CREATE_INDEX=true \\
		INPUT=$inputDir/\$INPUT_FILE OUTPUT=$indexedBamDir/\$OUTPUT_FILE
	EOF
	$indexBams =~ s/\t//g;
	push @jobScripts, $indexBams;
	$inputDir = $indexedBamDir;
}

# Step 5: Indel Realignment for each BAM file
if ($stepsToRun{'5'}) {
	unless ($stepsToRun{'4'}) {
		my @bamFileCheck = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
		my @indexFileCheck = split /\n/, `ls $inputDir/*bai 2> /dev/null`;
		$numFiles = scalar @bamFileCheck;
		unless ($numFiles > 0) {
			print "Error: No files ending in \".bam\" found in directory $inputDir.\n";
			exit;
		}
		unless (scalar @indexFileCheck == $numFiles) {
			print "Error: It appears that there are not an equal number of BAM files and index files (.bai) in $inputDir.\n";
			exit;
		}
	}

	my $realignDir = $outDir . '/realigned_bam';
	my $logDir = $realignDir . '/logs';
	push @directoriesToMake, ($realignDir, $logDir);

	# indel realignment involves running GATK's RealignerTargetCreator and IndelRealigner
	my $indelRealign = << "	EOF";
		\#\$ -N realign -j y -o $logDir/realign.\$TASK_ID.txt -V -S /bin/bash -t 1-$numFiles -tc $maxJobs -l mem_free=8G -pe parallel 4
		cd $inputDir
		BAM_FILES=(*bam)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${BAM_FILES[\$INDEX]}
		THIS_LOG=\${INPUT_FILE/%_indexed.bam/.txt}
		THIS_LOG=\${THIS_LOG/%bam/txt} # in case input file doesn't have "_indexed" in it
		mv "$logDir/realign.\$SGE_TASK_ID.txt" "$logDir/realign.\$SGE_TASK_ID.\$THIS_LOG"
		OUTPUT_FILE=`echo \$INPUT_FILE | sed -e "s/\\.bam/_realigned.bam/"`
		JOB_NUM=\$((SGE_TASK_ID))
		WORK_DIR="$realignDir/work/\$JOB_NUM"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		cd \$WORK_DIR
		echo "Creating target interval list for \$INPUT_FILE"
		##REPORT
		#Indel Realignment:
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk  -T RealignerTargetCreator -nt 4 $intervals $padding -isr INTERSECTION -R $refGenome \\
		-I $inputDir/\$INPUT_FILE -known $goldenIndels -o \$WORK_DIR/target_intervals.list --filter_mismatching_base_and_quals
		##END_REPORT
		echo "Now doing the actual realignment"
		##REPORT
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T IndelRealigner -R $refGenome -I $inputDir/\$INPUT_FILE \\
		-targetIntervals \$WORK_DIR/target_intervals.list -known $goldenIndels -o $realignDir/\$OUTPUT_FILE --filter_mismatching_base_and_quals
		##END_REPORT
	EOF

	$indelRealign =~ s/\t//g;
	push @jobScripts, $indelRealign;
	$inputDir = $realignDir;

}

# Step 6: Base quality score recalibration
if ($stepsToRun{'6'}) {
	unless ($stepsToRun{'5'}) {
		my @bamFileCheck = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
		my @indexFileCheck = split /\n/, `ls $inputDir/*bai 2> /dev/null`;
		$numFiles = scalar @bamFileCheck;
		unless ($numFiles > 0) {
			print "Error: No files ending in \".bam\" found in directory $inputDir.\n";
			exit;
		}
		unless (scalar @indexFileCheck == $numFiles) {
			print "Error: It appears that there are not an equal number of BAM files and index files (.bai) in $inputDir.\n";
			exit;
		}
	}

	my $recalDir = $outDir . '/recal_bam';
	my $logDir = $recalDir . '/logs';
	push @directoriesToMake, ($recalDir, $logDir);

	my $baseRecal = << "	EOF";
		\#\$ -N bqRecal -j y -o $logDir/bqsr.\$TASK_ID.txt -V -S /bin/bash -t 1-$numFiles -tc $maxJobs -l mem_free=8G -pe parallel 4
		cd $inputDir
		BAM_FILES=(*bam)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${BAM_FILES[\$INDEX]}
		THIS_LOG=\${INPUT_FILE/%_indexed_realigned.bam/.txt}
		THIS_LOG=\${THIS_LOG/%bam/txt}
		mv "$logDir/bqsr.\$SGE_TASK_ID.txt" "$logDir/bqsr.\$SGE_TASK_ID.\$THIS_LOG"
		OUTPUT_FILE=`echo \$INPUT_FILE | sed -e "s/\\.bam/_recal.bam/"`
		JOB_NUM=\$((SGE_TASK_ID))
		WORK_DIR="$recalDir/work/\$JOB_NUM"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		cd \$WORK_DIR
		echo "Performing base quality score recalibration for \$INPUT_FILE"
		##REPORT
		#Base quality score recalibration:
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T BaseRecalibrator $intervals $padding -isr INTERSECTION -nct 4 \\
		-R $refGenome -I $inputDir/\$INPUT_FILE -knownSites $goldenIndels -knownSites $dbSNP138 -o \$WORK_DIR/recal_data.table
		##END_REPORT
		echo "Applying recalibration..."
		##REPORT
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T PrintReads -nct 4 -R $refGenome -I $inputDir/\$INPUT_FILE \\
		-BQSR \$WORK_DIR/recal_data.table -o $recalDir/\$OUTPUT_FILE
	EOF

	## Currently unused code for doing a second pase of BaseRecalibrator followed by before-and-after plots (takes too much time and disk space to do it again)
		# echo "Doing a second pass..."
		# 	time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T BaseRecalibrator $intervals $padding -isr INTERSECTION -nct 4 \\
		# 	-R $refGenome -I $inputDir/\$INPUT_FILE -knownSites $goldenIndels -knownSites $dbSNP138 -BQSR \$WORK_DIR/recal_data.table -o \$WORK_DIR/post_recal_data.table
		# 	echo "Making before-and-after plot..."
		# 	time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T AnalyzeCovariates -R $refGenome -before \$WORK_DIR/recal_data.table \\
		# 	-after \$WORK_DIR/post_recal_data.table -plots \$WORK_DIR/recalibration_plots.pdf

	$baseRecal =~ s/\t//g;
	push @jobScripts, $baseRecal;
	$inputDir = $recalDir;
}

# Step 7: Deal with samples with more than one BAM file
if ($stepsToRun{'7'}) {
	unless ($stepsToRun{'6'}) {
		my @bamFileCheck = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
		my @indexFileCheck = split /\n/, `ls $inputDir/*bai 2> /dev/null`;
		$numFiles = scalar @bamFileCheck;
		unless ($numFiles > 0) {
			print "Error: No files ending in \".bam\" found in directory $inputDir.\n";
			exit;
		}
		unless (scalar @indexFileCheck == $numFiles) {
			print "Error: It appears that there are not an equal number of BAM files and index files (.bai) in $inputDir.\n";
			exit;
		}
	}

	my $finalBamDir = $outDir . '/final_bam';
	my $logDir = $finalBamDir . '/logs';
	push @directoriesToMake, ($finalBamDir, $logDir);

	# Note: calling a subroutine to do the job submission work in real-time
	my $mergeBams = << "	EOF";
		\#\$ -N mergeSamples -j y -o $logDir/merge_control.log -V -S /bin/bash -l h=ihg-node-1|ihg-node-27|ihg-node-50
		##REPORT
		#Merging same-sample BAM files (if applicable):
		perl $0 -mergeSameSampleBams -in $inputDir -out $finalBamDir -max $maxJobs
		##END_REPORT
	EOF

	$mergeBams =~ s/\t//g;
	push @jobScripts, $mergeBams;
	$inputDir = $finalBamDir;
}

# Step 8: Run HaplotypeCaller to produce GVCF files
if ($stepsToRun{'8'}) {
	# Currently no pedigree file needed since the special X/Y processing isn't happening
	# unless ($pedFile) {
	# 	print "Error: For haplotype calling, you must supply a pedigree file using the '-ped' option.\n";
	# 	exit;
	# } else {
	# 	$pedFile = abs_path $pedFile or die $!;
	# }
	unless ($stepsToRun{'7'}) {
		my @bamFileCheck = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
		my @indexFileCheck = split /\n/, `ls $inputDir/*bai 2> /dev/null`;
		$numFiles = scalar @bamFileCheck;
		unless ($numFiles > 0) {
			print "Error: No files ending in \".bam\" found in directory $inputDir.\n";
			exit;
		}
		unless (scalar @indexFileCheck == $numFiles) {
			print "Error: It appears that there are not an equal number of BAM files and index files (.bai) in $inputDir.\n";
			exit;
		}
	}

	my $vcfDir = $outDir . '/vcf';
	my $logDir = $vcfDir . '/logs';
	push @directoriesToMake, ($vcfDir, $logDir);

	my $makeVCFs = << "	EOF";
		\#\$ -N makeVCF -j y -o $logDir/hc.\$TASK_ID.txt -V -S /bin/bash -t 1-$numFiles -tc $maxJobs -l mem_free=8G
		cd $inputDir
		BAM_FILES=(*bam)
		INDEX=\$((SGE_TASK_ID-1))
		JOB_NUM=\$((SGE_TASK_ID))
		INPUT_FILE=\${BAM_FILES[\$INDEX]}
		OUTPUT_FILE=`echo \$INPUT_FILE | sed -e "s/\\.bam/.vcf/"`
		THIS_LOG=\${INPUT_FILE/%bam/txt}
		mv "$logDir/hc.\$SGE_TASK_ID.txt" "$logDir/hc.\$SGE_TASK_ID.\$THIS_LOG"
		WORK_DIR="$vcfDir/work/\$JOB_NUM"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		cd \$WORK_DIR
		##REPORT
		#Haplotype calling:
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T HaplotypeCaller -R $refGenome -I $inputDir/\$INPUT_FILE \\
		-o $vcfDir/\$OUTPUT_FILE $intervals $padding -isr INTERSECTION -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --read_filter BadCigar \\
		--annotation StrandOddsRatio --annotation AlleleBalanceBySample --annotation DepthPerSampleHC --annotation MappingQualityZeroBySample \\
		--annotation StrandBiasBySample --annotation GenotypeSummaries
	EOF
	$makeVCFs =~ s/\t//g;
	push @jobScripts, $makeVCFs;
	$inputDir = $vcfDir;
}


# Step 9: CombineGVCFs
if ($stepsToRun{'9'}) {
	my $combineDir = $outDir . '/comb_gvcf';
	push @directoriesToMake, $combineDir;

	# Note that a subroutine is being called to submit the actual CombineGVCF jobs
	my $comb_master = << "	EOF";
		\#\$ -N comb_master -j y -o $combineDir/combine_master.log -V -S /bin/bash -l h=ihg-node-1|ihg-node-27|ihg-node-50
		cd $inputDir
		VCF_FILES=(*vcf)
		##REPORT
		#Combining VCFs into batches of 50 and running CombineGVCFs:
		perl $0 -in $inputDir -out $combineDir -combineGVCFs
	EOF
	$comb_master =~ s/\t//g;
	push @jobScripts, $comb_master;
	$inputDir = $combineDir;
}

# Step 10: Joint genotyping
if ($stepsToRun{'10'}) {
	my $pedigreeArgument = '';
	unless ($pedFile) {
		print "Warning: a pedigree file hasn't been specified, so there will be no inbreeding coefficient annotation in the joint VCF.\n";
	} else {
		$pedFile = abs_path $pedFile or die $!;
		$pedigreeArgument = "--pedigree $pedFile --annotation InbreedingCoeff";
	}

	unless ($stepsToRun{'9'} || $stepsToRun{'8'}) {
		my @vcfCheck = split /\n/, `ls $inputDir/*vcf 2> /dev/null`;
		$numFiles = scalar @vcfCheck;
		unless ($numFiles > 0) {
			print "Error: No files ending in \".vcf\" found in directory $inputDir.\n";
			exit;
		}
	}

	# adding an interval specification for when haplotype calling hasn't already included the interval
	my $intervalOptions = '';
	unless ($stepsToRun{'8'}) {
		$intervalOptions = "$intervals $padding -isr INTERSECTION";
	}

	my $jointDir = $outDir . '/joint';
	my $workDir = $jointDir . '/work'; # gets used in CombineGVCFs scenario
	my $logDir = $jointDir . '/logs';
	push @directoriesToMake, ($jointDir, $workDir, $logDir);

	my $fileEndingString = '(' . (join " ", @vcfEndings) . ')'; # see top of program

	# Note: Memory/threads may need to be altered for some jobs (larger ones need more memory, fewer threads; smaller jobs can run faster with more threads)
	my $jointGeno = << "	EOF";
		\#\$ -N jointGeno -j y -o $logDir/joint-\$TASK_ID.log -V -S /bin/bash -l mem_free=12G -pe parallel 6 -t 1-$numJointJobs -tc $maxJobs
		cd $inputDir
		FILE_GROUPS=$fileEndingString
		INDEX=\$((SGE_TASK_ID-1))
		CURR_GROUP=\${FILE_GROUPS[\$INDEX]}
		mv $logDir/joint-\$SGE_TASK_ID.log $logDir/joint-\$SGE_TASK_ID-\$CURR_GROUP.log
		WORK_DIR="$workDir/\$SGE_TASK_ID"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp
		VCF_FILES=(*\$CURR_GROUP)
		VAR_LINE=''
		for FILE in "\${VCF_FILES[@]}"
		do
			VAR_LINE+=" --variant $inputDir/\$FILE"
		done
		##REPORT
		#Joint genotyping:
		time $javaLoc/java -XX:+UseSerialGC -Xmx12g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T GenotypeGVCFs -nt 6 -R $refGenome \$VAR_LINE \\
		--out $jointDir/joint-\$CURR_GROUP $intervalOptions $pedigreeArgument --annotation StrandOddsRatio --annotation BaseQualityRankSumTest \\
		--annotation ChromosomeCounts --annotation Coverage --annotation FisherStrand \\
		--annotation MappingQualityRankSumTest --annotation MappingQualityZero \\
		--annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation VariantType \\
		--annotation DepthPerAlleleBySample --annotation AlleleBalanceBySample --annotation MappingQualityZeroBySample \\
		--annotation StrandBiasBySample --annotation DepthPerSampleHC --annotation GenotypeSummaries \\
	EOF
	$jointGeno =~ s/\t//g;
	push @jobScripts, $jointGeno;
	$inputDir = $jointDir;
}

# Step 11: SNP recalibration
if ($stepsToRun{'11'}) {
	unless ($stepsToRun{'10'}) {
		my @jointCheck = split /\n/, `ls $inputDir/*vcf 2> /dev/null`;
		$numJointJobs = scalar @jointCheck; # can overwrite the expected number of jobs if a different number of files is found
		unless ($numJointJobs > 0) {
			print "Error: There should be at least one joint VCF file ending in \".vcf\" in your input folder.\n";
			exit;
		}
	}

	my $snpRecalDir = $outDir . '/snp_recal';
	my $logDir = $snpRecalDir . '/logs';
	my $workDir = $snpRecalDir . '/work';
	push @directoriesToMake, ($snpRecalDir, $logDir, $workDir);

	# Note: Matching GATK's argument recommendations for recalibration on whole genome variants (just adding DP annotation)
	my $snpRecal = << "	EOF";
		\#\$ -N snp_recal -j y -o $logDir/snpRecal-\$TASK_ID.log -V -S /bin/bash -l mem_free=16G -t 1-$numJointJobs -tc $maxJobs
		cd $inputDir
		ALL_FILES=(*vcf)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${ALL_FILES[\$INDEX]}
		TMP_DIR=$workDir/\$SGE_TASK_ID/tmp
		mkdir -p \$TMP_DIR
		OUTPUT_FILE=\${INPUT_FILE/.vcf/_snp-recal.vcf}
		BASENAME=\${INPUT_FILE/.vcf/}
		mv $logDir/\$SGE_TASK_ID-snpRecal.log $logDir/snpRecal-\$SGE_TASK_ID-\$BASENAME.log
		##REPORT
		SNP recalibration:
		time $javaLoc/java -XX:+UseSerialGC -Xmx16g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T VariantRecalibrator -R $refGenome -input \$INPUT_FILE \\
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $recalTraining1 \\
		-resource:omni,known=false,training=true,truth=true,prior=12.0 $recalTraining2 \\
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 $recalTraining3 \\
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $recalTraining4 \\
		-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff -mode SNP \\
		-recalFile $snpRecalDir/\$BASENAME-SNP.recal -tranchesFile $snpRecalDir/\$BASENAME-SNP.tranches -rscriptFile $snpRecalDir/\$BASENAME-SNP_recal_plots.R
		##END_REPORT
		echo "Now applying recalibration..."
		##REPORT
		time $javaLoc/java -XX:+UseSerialGC -Xmx16g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T ApplyRecalibration -R $refGenome -input \$INPUT_FILE \\
		-recalFile $snpRecalDir/\$BASENAME-SNP.recal -tranchesFile $snpRecalDir/\$BASENAME-SNP.tranches -o $snpRecalDir/\$OUTPUT_FILE \\
		 -mode SNP --ts_filter_level 99.5
	EOF
	$snpRecal =~ s/\t//g;
	push @jobScripts, $snpRecal;
	$inputDir = $snpRecalDir;
}

# Step 12: Indel recalibration
if ($stepsToRun{'12'}) {
	unless ($stepsToRun{'11'}) {
		my @jointCheck = split /\n/, `ls $inputDir/*vcf 2> /dev/null`;
		$numJointJobs = scalar @jointCheck; # can overwrite the expected number of jobs if a different number of files is found
		unless ($numJointJobs > 0) {
			print "Error: There should be at least one joint VCF file ending in \".vcf\" in your input folder.\n";
			exit;
		}
	}

	my $indelRecalDir = $outDir . '/indel_recal';
	my $logDir = $indelRecalDir . '/logs';
	my $workDir = $indelRecalDir . '/work';
	push @directoriesToMake, ($indelRecalDir, $logDir);

	# Note: Matching GATK's argument recommendations for recalibration on whole genome variants (just adding DP annotation)
	my $indelRecal = << "	EOF";
		\#\$ -N indel_recal -j y -o $logDir/indelRecal-\$TASK_ID.log -V -S /bin/bash -l mem_free=16G -t 1-$numJointJobs -tc $maxJobs
		cd $inputDir
		ALL_FILES=(*vcf)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${ALL_FILES[\$INDEX]}
		TMP_DIR=$workDir/\$SGE_TASK_ID/tmp
		mkdir -p \$TMP_DIR
		OUTPUT_FILE=\${INPUT_FILE/.vcf/_indel-recal.vcf}
		BASENAME=\${INPUT_FILE/_snp-recal.vcf/}
		BASENAME=\${BASENAME/.vcf/}
		mv $logDir/indelRecal-\$SGE_TASK_ID.log $logDir/indelRecal-\$SGE_TASK_ID-\$BASENAME.log
		##REPORT
		#Indel recalibration:
		time $javaLoc/java -XX:+UseSerialGC -Xmx16g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T VariantRecalibrator -R $refGenome -input \$INPUT_FILE --maxGaussians 4 \\
		-resource:mills,known=false,training=true,truth=true,prior=12.0 $goldenIndels \\
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP138 \\
		-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -mode INDEL \\
		-recalFile $indelRecalDir/\$BASENAME-INDEL.recal -tranchesFile $indelRecalDir/\$BASENAME-INDEL.tranches -rscriptFile $indelRecalDir/\$BASENAME-INDEL-recal_plots.R
		##END_REPORT
		echo "Now applying recalibration..."
		##REPORT
		time $javaLoc/java -XX:+UseSerialGC -Xmx16g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T ApplyRecalibration -R $refGenome -input \$INPUT_FILE \\
		-recalFile $indelRecalDir/\$BASENAME-INDEL.recal -tranchesFile $indelRecalDir/\$BASENAME-INDEL.tranches -o $indelRecalDir/\$OUTPUT_FILE \\
		 -mode INDEL --ts_filter_level 90.0
	EOF
	$indelRecal =~ s/\t//g;
	push @jobScripts, $indelRecal;
	$inputDir = $indelRecalDir;
}

# Step 13: De novo variant calling (BETA)
if ($stepsToRun{'13'}) {
	unless ($pedFile) {
		print "Error: For de novo calling, you must supply a pedigree file using the '-ped' option.\n";
		exit;
	} else {
		$pedFile = abs_path $pedFile or die $!;
	}

	unless ($phenFile) {
		print "Error: For de novo calling, you must supply a phenotype file using the '-phe' option.\n";
		print "See the PlinkSeq website for information on proper formatting.\n";
		exit;
	} else {
		$phenFile = abs_path $phenFile or die $!;
	}

	unless ($stepsToRun{'12'}) {
		my @jointCheck = split /\n/, `ls $inputDir/*vcf 2> /dev/null`;
		$numJointJobs = scalar @jointCheck; # can overwrite the expected number of jobs if a different number of files is found
		unless ($numJointJobs > 0) {
			print "Error: For de novo variant calling, there should at least one VCF file ending in \".vcf\" in your input folder.\n";
			exit;
		}
	}

	my $projectName;
	print "Enter a name for your PlinkSeq project: ";
	while (my $userInput = <STDIN>) {
		chomp $userInput;
		$projectName = $userInput;
		last if ($projectName =~ /[\w-]+/);
		print "Please enter a more normal-seeming project name: ";
	}

	my $deNovoDir = $outDir . '/dn_calling';
	my $logDir = $deNovoDir . '/logs';
	my $scratchDir = $deNovoDir . '/work';
	push @directoriesToMake, ($deNovoDir, $logDir, $scratchDir);

	# The awk line generates a PlinkSeq format pedigree file from a regular pedigree file
	my $callDeNovo = << "	EOF";
		\#\$ -N dn_call -j y -o $logDir/pseq-\$TASK_ID.log -V -S /bin/bash -l mem_free=16G -t 1-$numJointJobs -tc $maxJobs
		cd $inputDir
		ALL_FILES=(*.vcf)
		INDEX=\$((SGE_TASK_ID-1))
		INPUT_FILE=\${ALL_FILES[\$INDEX]}
		cd $deNovoDir
		SCRATCH=$scratchDir/\$SGE_TASK_ID
		mkdir -p \$SCRATCH
		PS_PED=`basename $pedFile`
		PS_PED=$deNovoDir/converted_\$PS_PED
		if [ \$SGE_TASK_ID -eq "1" ]
		then
			# technically, there could be a problem if converted file isn't generated before another job in the array calls it, but seems unlikely
			awk -F "\\t" 'BEGIN {OFS="\\t"} {print \$2,\$1,\$2,\$3,\$4,\$5}' $pedFile >  \$PS_PED # PlinkSeq requires a special pedigree format!
		fi
		OUTPUT=\${INPUT_FILE/.vcf/}
		OUTPUT=\${OUTPUT/_indel-recal/}
		OUTPUT=\${OUTPUT/_snp-recal/}
		PROJECT=$projectName\_\$OUTPUT
		mv $logDir/pseq-\$SGE_TASK_ID.log $logDir/pseq-\$SGE_TASK_ID-\$OUTPUT.log
		##REPORT
		#PlinkSeq:
		$plinkSeq \$PROJECT new-project --resources $pseqGenome --scratch \$SCRATCH --vcf $inputDir/\$INPUT_FILE
		$plinkSeq \$PROJECT load-vcf
		TAG=\$PROJECT
		TAG+=_dn
		$plinkSeq \$PROJECT tag-file --id 1 --name \$TAG
		$plinkSeq \$PROJECT var-summary
		$plinkSeq \$PROJECT load-pedigree --file \$PS_PED
		$plinkSeq \$PROJECT load-pheno --file $phenFile
		$plinkSeq \$PROJECT denovo --mask any.filter.ex aac=1-5 --minHet_AB_alt 0.3 --minHet_AB_ref 0.3 --minHomRef_AB_ref 0.95 \\
		--minChildPL 20 --minParPL 20 --minMQ 30 --minChildDP 10 --minParDP 10 --out $deNovoDir/\$PROJECT-DP10.scan
		$plinkSeq \$PROJECT denovo --printTransmission --mask any.filter.ex --out  $deNovoDir/\$PROJECT-withTransmission_noFilters.scan
		$plinkSeq \$PROJECT i-stats --mask any.filter.ex
	EOF
	$callDeNovo =~ s/\t//g;
	push @jobScripts, $callDeNovo;
	$inputDir = $deNovoDir;
}

checkForExistingDirs(@directoriesToMake); # dealing with any directories from previous pipeline runs
for (@directoriesToMake) {
	unless (mkdir $_) {
		print "Failed to create directory $_.\n";
		print "Exiting. No jobs were submitted.\n";
		exit;
	}
}

# creating (or re-opening) a file that records run-time arguments
my $runSummary = "$outDir/pipe_run_summary.txt";
$inputDir = $initialInputDir; # resetting this so it prints correctly in log (clumsy, I know)
my $addSpace = (-e $runSummary) ? 1 : 0;
open my $LOG, '>>', $runSummary or die $!;
print $LOG "\n\n" if ($addSpace);
my @timeArray = localtime;

my ($hours, $minutes, $seconds) = @timeArray[2, 1, 0];
for ($minutes, $seconds) {
	$_ = "0$_" if ($_ < 10); # so 11:02:09 is not rendered as 11:2:9 (lame, right?)
}
my $time = join ":", ($hours, $minutes, $seconds);
$timeArray[4] += 1; # to get a normal month (default is 0=January, 11=December)
$timeArray[5] += 1900; # to get a normal year (default is years since 1900)
my $date = join "/", @timeArray[4, 3, 5];
print $LOG "Program: $0\nDate: $date, $time\n";
print $LOG "Working directory: " . getcwd() . "\n";
print $LOG "Command-line arguments:\n";
for (sort { length $a <=> length $b} keys %options) {
	my $bufferLength = 20;
	my $value = ${$options{$_}};
	next unless ($value);
	(my $option = $_) =~ s/[:|].*//;
	$bufferLength -= length $option;
	my $spacing = ' ' x $bufferLength;
	print $LOG "\t", $option, $spacing, $value, "\n";
}
print $LOG "\n";

# Submitting all the jobs (with holds so that they run one step at a time)
my ($previousJobID, $holdJob) = ('', '');
for (@jobScripts) {
	chomp ($previousJobID = `qsub -terse $holdJob <<"EOF"\n$_\nEOF`);
	print $LOG getLinesForReport($_), "\n\n";
	if (${^CHILD_ERROR_NATIVE}) {
		print "Warning: At least one job didn't seem to submit correctly, and any subsequent jobs were not submitted.\nCheck with \"qsub\" to see what's running.\n";
		print $LOG "---An error occurred here---\n";
		print $LOG "\n", ('-') x 50, "\n";
		exit;
	}
	$previousJobID =~ s/\..*//; # for jobs with array job suffixes, removes the suffix to get just the job number
	$holdJob = "-hold_jid $previousJobID";
}
print $LOG "\n", ('-') x 50, "\n";
close $LOG;

print "All jobs have have been submitted to the cluster. Check their progress with \"qstat.\"\n";
print "When the jobs are finished, look for output files in your specified output directory.\n";




## Subroutines ##

# This subroutine reads a fastq file along with the Fastq-SampleID list file to construct a header line for input into BWA
sub getHeaderFromFq($$) {
	my ($fq, $sampleFile) = @_;
	(my $basename = $fq) =~ s/.*\///;
	$basename =~ s/\.(fastq|fq)(\.gz)?$//; # getting rid of extensions (easier to compare with sample list)
	# sample file should already have been checked for errors, so not checking now
	my $sampleID;
	open my $IN, '<', $sampleFile or die $!;
	while (my $line = <$IN>) {
		chomp $line;
		next unless ($line);
		my @fields = split /\t/, $line;
		my $file = $fields[0];
		$file =~ s/.*\///; # strip path
		$file =~ s/\.(fastq|fq)(\.gz)?$//; # getting rid of file extensions (and on the next line, split fastq numbering)
		$file =~ s/(_\d{3}$)// unless ($file eq $basename); # the "unless" is for the situation where the SampleID ends in an underscore followed by three digits (just like the split fq numbering)
		if ($file eq $basename) {
			$sampleID = $fields[1];
			last;
		}
	}
	close $IN;
	# Signaling an error when no Sample ID is found (likely because the sample is missing from the Fastq-SampleID list)
	unless ($sampleID) {
		print 'ERROR';
		exit;
	}
	chomp (my $headerLine = `less -f $fq | head -1`); # less is needed because file may be gzipped
	$headerLine =~ s/^@//; # gets in the way
	my @fields = split /:/, $headerLine;
	my $ID = 'ID:' . $fields[2] . ':' . $fields[3];
	my $SM = 'SM:' . $sampleID;
	my $PL = 'PL:ILLUMINA';
	my $LB = 'LB:lib.' . $sampleID . '.' . $fields[2];
	my $PU = 'PU:' . (join ".", ($fields[0], $fields[2], $fields[1], $fields[3], $fields[-1]));
	my $readGroupInfo = "'" . (join '\t', ('@RG', $ID, $SM, $PL, $LB, $PU)) . "'";
	print $readGroupInfo;
}

# This subroutine orchestrates the task of merging same-sample BAM files (and repeating indel realignment and MarkDuplicates)
sub mergeBamsAndRealign($$) {
	my ($previousJobID, $holdJob) = ('', '');
	my ($inputDir, $outDir) = @_;
	$inputDir = abs_path $inputDir;
	$outDir = abs_path $outDir;
	my $logDir = $outDir . '/logs';
	mkdir $logDir; # should already exist

	# getting SampleIDs of input BAM files
	my @bams = split /\n/, `ls $inputDir/*bam 2> /dev/null`;
	my %sampleIDsToFiles;
	for (@bams) {
		my $sampleID = getSampleIDFromBam($_);
		push @{$sampleIDsToFiles{$sampleID}}, $_;
	}

	# constructing gatk commands to merge groups of BAM files that share a SampleID
	my @mergeCommands;
	my @outputNames;
	for (keys %sampleIDsToFiles) {
		# when just one file corresponds to a SampleID, adding a link to the file in the final BAM directory
		if (scalar @{$sampleIDsToFiles{$_}} == 1) {
			my $file = shift @{$sampleIDsToFiles{$_}};
			my $basename = $_ . '.bam'; # provisionally using Sample ID to name BAM files
			(my $indexFile = $file) =~ s/\.bam$/.bai/;
			# creating hard links instead of symbolic ones so that the recalibrated BAM folder can be safely deleted later
			my $indexBasename = $_ . '.bai';
			`ln -T $file $outDir/$basename`;
			`ln -T $indexFile $outDir/$indexBasename`; # got to link the index file too
		} else {
			my $mergeCommand;
			for my $bam (@{$sampleIDsToFiles{$_}}) {
				$mergeCommand .= "-I $bam ";
			}
			my $mergedBam = "$_.bam"; # just using the SampleID for the file name
			push @outputNames, $mergedBam;
			print "Files to be merged for sample $_: $mergeCommand\n";
			push @mergeCommands, $mergeCommand;                                                        
		}
	}
	exit unless (@mergeCommands); # if nothing to merge, done with this instance of the script

	# need a file listing the commands in order to handle very large jobs
	my $commandFile = "$logDir/mergeCommands.txt";
	open my $OUT, '>', $commandFile or die $!;
	print $OUT $_, "\n" for (@mergeCommands);
	close $OUT;

	my $outputFileNames = "$logDir/outputFileNames.txt";
	open $OUT, '>', $outputFileNames or die $!;
	print $OUT $_, "\n" for (@outputNames);
	close $OUT;

	my $numMerges = scalar @mergeCommands;
	my $mergeBams = <<"	EOF";
		\#\$ -N mergeBam -j y -o $logDir/realign.\$TASK_ID.txt -V -S /bin/bash -t 1-$numMerges -tc $maxJobs -l mem_free=8G
		INDEX=\$((SGE_TASK_ID))
		CURR_COMMAND=`awk "NR==\$INDEX" $commandFile`
		CURR_OUTPUT=`awk "NR==\$INDEX" $outputFileNames`
		THIS_LOG=\${CURR_OUTPUT/%bam/txt}
		mv "$logDir/realign.\$SGE_TASK_ID.txt" "$logDir/realign.\$SGE_TASK_ID.\$THIS_LOG"
		MERGE_OUT=`echo \$CURR_OUTPUT | sed -e "s/\\.bam/_initial.bam/"`
		DEDUP_OUT=`echo \$CURR_OUTPUT | sed -e "s/\\.bam/_dedup.bam/"`
		METRICS_FILE=`echo \$CURR_OUTPUT | sed -e "s/\\.bam/_metrics.txt/"`
		JOB_NUM=\$((SGE_TASK_ID))
		WORK_DIR="$outDir/work/\$JOB_NUM"
		PROC_DIR="$outDir/work/"
		mkdir -p \$WORK_DIR
		TMP_DIR=\$WORK_DIR/tmp # to avoid overwriting files across jobs
		cd \$WORK_DIR
		echo \$CURR_COMMAND
		FULL_COMMAND="time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T PrintReads -R $refGenome \$CURR_COMMAND -o \$PROC_DIR/\$MERGE_OUT"
		echo \$FULL_COMMAND
		eval \$FULL_COMMAND
		echo "Marking duplicates..."
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $picardLoc/MarkDuplicates.jar METRICS_FILE=\$PROC_DIR/\$METRICS_FILE CREATE_INDEX=true \\
		INPUT=\$PROC_DIR/\$MERGE_OUT OUTPUT=\$PROC_DIR/\$DEDUP_OUT
		rm \$PROC_DIR/\$MERGE_OUT # JM: REMOVE THIS LATER
		echo "Generated \$PROC_DIR/\$DEDUP_OUT and deleted \$PROC_DIR/\$MERGE_OUT"
		echo "Creating target interval list for \$PROC_DIR/\$DEDUP_OUT"
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk  -T RealignerTargetCreator $intervals $padding -isr INTERSECTION \\
		-R $refGenome -I \$PROC_DIR/\$DEDUP_OUT -known $goldenIndels -o \$WORK_DIR/target_intervals.list
		echo "Now doing the actual indel realignment..."
		time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T IndelRealigner -R $refGenome -I \$PROC_DIR/\$DEDUP_OUT \\
		-targetIntervals \$WORK_DIR/target_intervals.list -known $goldenIndels -o $outDir/\$CURR_OUTPUT
	EOF

	$mergeBams =~ s/\t//g;
	chomp ($previousJobID = `qsub -terse $holdJob <<"EOF"\n$mergeBams\nEOF`);
	$previousJobID =~ s/\..*//; # removing the array job info to get just the job number
	$holdJob = "-hold_jid $previousJobID";

	my $markDone = <<"	EOF";
		\#\$ -N marker -j y -o /dev/null -V -S /bin/bash
		touch $outDir/.done
	EOF

	$markDone =~ s/\t//g;
	`qsub -terse $holdJob <<"EOF"\n$markDone\nEOF`;

	# going to tread water until the merging is done so that the pipeline doesn't continue without it
	while () {
		last if (-e "$outDir/.done");
		sleep 15;
	}
	
}

# This subroutine divides a large group of VCFs into batches of 50 and sumbits a CombineGVCF job for each batch
sub combineGVCFs($$) {
	my ($previousJobID, $holdJob) = ('', '');
	my ($inputDir, $outDir) = @_;
	my @vcfs = split /\n/, `ls *vcf 2> /dev/null`;
	my %fileGroups;
	my $numFileGroups = 0;
	my $maxVCF = 1; # changed from 50 to 90 to 25
	while (@vcfs) {
		$numFileGroups++ if ($numFileGroups == 0 || scalar @{$fileGroups{$numFileGroups}} >= $maxVCF);
		my $currentVCF = shift @vcfs;
		push @{$fileGroups{$numFileGroups}}, $currentVCF;
	}

	my $commandList = $outDir . '/file_batch_list.txt';
	my $logDir = $outDir . '/logs';
	mkdir $logDir;

	open my $OUT, '>', $commandList or die $!;
	# going to print out entire groups of VCFs on each line of list file
	for my $currentGroup (sort {$a <=> $b} keys %fileGroups) {
		for my $file (@{$fileGroups{$currentGroup}}) {
			(my $basename = $file) =~ s/.*\///;
			print $OUT "--variant $basename ";
		}
		print $OUT "\n";
	}

	my $chrGroupsList = '(' . (join ' ', @chrGroups) . ')'; # see top of program
	my $numJobs = $numFileGroups * $numChrGroups;

	my $combineGVCFs = <<"	EOF";
		\#\$ -N comb_gvcf -j y -o $logDir/combine.\$TASK_ID.txt -V -S /bin/bash -t 1-$numJobs -tc $maxJobs -l mem_free=8G
		cd $inputDir
		INDEX=\$((SGE_TASK_ID))
		CHR_GROUPS=$chrGroupsList
		CHR_INDEX=\$(((\$INDEX-1)\%$numChrGroups))
		CURR_INTERVAL=\${CHR_GROUPS[\$CHR_INDEX]}
		FILE_INDEX=\$((((\$INDEX-1)/$numChrGroups)+1))
		VARIANTS=`awk "NR==\$FILE_INDEX" $commandList`
		INTERVAL_NAME_TEMP=\${CURR_INTERVAL//-L /}
		INTERVAL_NAME=\${INTERVAL_NAME_TEMP// /,}
		OUTPUT="$outDir/batch\$FILE_INDEX-chr\$INTERVAL_NAME.vcf"
		mv $logDir/combine.\$SGE_TASK_ID.txt $logDir/comb-\$SGE_TASK_ID-batch\$FILE_INDEX-chr\$INTERVAL_NAME.txt
		TMP_DIR="$outDir/work/\$INDEX/tmp"
		mkdir -p \$TMP_DIR
		COMMAND="time $javaLoc/java -XX:+UseSerialGC -Xmx8g -Djava.io.tmpdir=\$TMP_DIR -jar $gatk -T CombineGVCFs \$CURR_INTERVAL -R $refGenome \$VARIANTS -o \$OUTPUT"
		echo -e "Running this:\n\$COMMAND"
		eval \$COMMAND
	EOF

	$combineGVCFs =~ s/\t//g;
	chomp ($previousJobID = `qsub -terse $holdJob <<"EOF"\n$combineGVCFs\nEOF`);
	$previousJobID =~ s/\..*//; # removing the array job info to get just the job number
	$holdJob = "-hold_jid $previousJobID";

	my $markDone = <<"	EOF";
		\#\$ -N marker -j y -o /dev/null -V -S /bin/bash
		touch $outDir/.done
	EOF

	$markDone =~ s/\t//g;
	`qsub -terse $holdJob <<"EOF"\n$markDone\nEOF`;

	# going to tread water until the merging is done so that the pipeline doesn't continue without it
	while () {
		last if (-f "$outDir/.done");
		sleep 15;
	}
}

# This subroutine uses Samtools to get the SampleID from a BAM file's header (Note: Might behave badly if a file has more than one Sample ID)
sub getSampleIDFromBam($) {
	my $bam = shift @_;
	my $readGroupLine = `$samtoolsLoc view -H $bam | grep -m 1 '^\@RG.*SM:'`;
	(my $sampleID = $readGroupLine) =~ s/.*\tSM:([^\t]*).*/$1/;
	chomp $sampleID;
	return $sampleID;
}

# This subroutine determines whether a directory with a given path exists, and tries to remove it if it does (with permission if directory is not empty)
sub checkForExistingDirs(@) {
	my @directories = @_;
	my @conflictingDirs;
	for my $dir (@directories) {
		if (-f $dir) {
			print "There is a file at $dir, but this pipeline needs to create a directory with that name. Please move or delete that file.\n";
			print "Exiting. No jobs were submitted.\n";
			exit;
		}
		# removes empty directory trees
		if (-d $dir) {
			finddepth(sub{rmdir}, $dir);
			rmdir $dir;
		}
		push @conflictingDirs, $dir if (-d $dir);
	}
	if (@conflictingDirs) {
		print "\nThe pipeline needs to create and use the following directories, but they already exist:\n";
		print $_, "\n" for (@conflictingDirs);
		print "Is it okay to overwrite the contents of these directories (yes/no)? ";
		while () {
			chomp (my $answer = <STDIN>);
			if ($answer =~ /^no$/i) {
				print "Exiting. No jobs were submitted.\n";
				exit;
			}
			elsif ($answer =~ /^yes$/i)  {
				for (@conflictingDirs) {
					`rm -R $_` if (-d $_); # to avoid getting an error messaging when attempting to remove a directory whose parent directory has already been removed
				}
				last;
			} else {
				print "Invalid response. Answer \"yes\" to delete the directories, or \"no\" to exit: ";
			}
		}
	}
}

# This subroutine parses heredoc bash scripts for the lines that should go into the run report file (if all lines were printed, it would be confusing)
sub getLinesForReport($) {
	my $script = shift @_;
	$script =~ s/\\\n//gs;
	my @lines = split /\n/, $script;
	my @reportLines;
	my $report = 0;
	for (@lines) {
		chomp $_;
		if ($_ eq '##REPORT') {
			$report = 1;
			next;
		}
		if ($_ eq '##END_REPORT') {
			$report = 0;
			next;
		}
		$_ =~ s/^#//;
		push @reportLines, $_ if ($report);
	}
	my $output = join "\n", @reportLines;
	return $output;
}
# This subroutine prints the steps in the pipeline and asks the user which ones to run
sub getPipelineChoice() {
	my $msg = <<"	EOF";

		This pipeline performs the following steps:
		1) Merge split Fastq files and standardize file extensions
		2) Run BWA on Fastq files to produce SAM files (files must end in .fastq.gz or .fastq)
		3) Convert SAM files to BAM files and sort BAM files by coordinate
		4) Mark duplicates in sorted BAM files and produce BAM index files
		5) Indel realignment
		6) Base quality score recalibration
		7) Merge same-sample BAMs and re-run steps 4-5 on merged BAMs
		8) Run HaplotypeCaller on sorted, indexed BAM files to produce gVCF files
		9) CombineGVCFs (optional, but highly recommended if you have >200 samples)
		10) GenotypeGVCFs (joint genotyping)
		11) SNP recalibration of joint VCF
		12) Indel recalibration of joint VCF
		13) De novo calling with PlinkSeq (beta)
		
	EOF

	$msg =~ s/\t//g;
	print $msg;
	my $numSteps = 13; # edit if pipeline changed
	my $combineStep = 9; # edit if pipeline changed
	print "Which step (1-$numSteps) do you wish to start with? ";
	my ($begin, $end);
	while ($begin = <STDIN>) {
		chomp $begin;
		last if ($begin =~ /^\d+$/ && $begin <= $numSteps);
		print "Invalid response. Enter the number from 1 to $numSteps that corresponds to the step to start with: ";
	}

	if ($begin == $numSteps) {
		$end = $numSteps;
	} else {
		print "Which step ($begin-$numSteps) do you wish to end with? ";
		while ($end = <STDIN>) {
			chomp $end;
			last if ($end =~ /^\d+$/ && $end >= $begin && $end <= $numSteps);
			print "Invalid response. Enter the number from $begin to $numSteps that corresponds to the step to end with: ";
		}
	}

	my @steps = (1 .. $numSteps);
	my %stepsToRun;
	$stepsToRun{$steps[$_ - 1]} = 1 for ($begin..$end);

	# special check for CombineGVCF
	if ($stepsToRun{$combineStep} && $stepsToRun{$combineStep - 1} && $stepsToRun{$combineStep + 1}) { # no need to ask if user explicitly specified step as first or last
		print "Would you like to combine VCFs before joint genotyping? This is optional, but strongly recommended " . 
			"if you have greater than 200 samples. Answer \"yes\" or \"no:\" ";
		while (my $answer = <STDIN>) {
			chomp $answer;
			if ($answer =~ /^y(es)?$/i) {
				print "Okay, CombineGVCFs will run.\n";
				last;
			} elsif ($answer =~ /^no?$/i) {
				delete $stepsToRun{$combineStep};
				print "Okay, CombineGVCFs will be skipped.\n";
				last;
			}
			else {
				print "Please answer \"yes\" or \"no\": ";
			}
		}
	}
	return %stepsToRun;
}

# This subroutine prints the documentation for this program
sub help() {
	my $message = <<"	EOF";
		
		Usage: perl $0 -in [dir] -out [dir]

		Description:
		This program processes fastq, SAM, BAM, and VCF files. It's possible to run the full Fastq -> GVCF pipeline,
		or to run any subset of steps in the process. Just run the program with appopriate arguments, and you'll
		be prompted to pick which type of processing to perform.

		Arguments:
		-in: specify the directory containing the files to be processed
		-out: give a name or path for a new output directory, or specify an existing one
		-list: specifies a file listing Sample IDs for each Fastq file (required for running BWA)
		-ped: specifies a pedigree file (required for PlinkSeq, recommended for joint genotyping)
		-interval: specifies an interval file for GATK (to specify multiple interval files, separate them by commas)
		-phenotype: specifies a phenotype file for PlinkSeq (only required for PlinkSeq)
		-maxjobs: the number of jobs to submit at once with qsub (default 91)
		-help: prints this message

		Pipeline Steps:
		1) Merge split Fastq files and standardize file extensions
		2) Run BWA on Fastq files to produce SAM files
		3) Convert SAM files to BAM files and sort BAM files by coordinate
		4) Mark duplicates in sorted BAM files and produce BAM index files
		5) Indel realignment
		6) Base quality score recalibration
		7) Deal with samples that have more than one BAM (see below)
		8) Run GATK's HaplotypeCaller, indexed BAM files to produce GVCF files
		9) CombineGVCFs(optional, but highly recommended when including >200 VCFs in analysis)
		10) GenotypeGVCFs (combine GVCFs into one joint VCF)
		11) SNP recalibration of joint VCF
		12) Indel recalibration of joint VCF
		13) PlinkSeq (beta)

		File Formats: Fastq files supplied as input to Step 1 of the pipeline must have filenames ending in .fq or .fastq,
		followed by .gz if they're compressed. Split  Fastq files (e.g. Samp1_R1_001.fq, Samp1_R1_002.fq, etc.) for any one
		sample should be either all compressed or all uncompressed. In addition, this pipeline assumes that all paired-end Fastq
		files have _R1 or _R2 preceding their file extensions, and that split Fastq files follow the format _R1_001, _R1_002,
		etc. Files without R1/R2 will be assumed to be interleaved. Interleaved and paired-end Fastq files cannot be included in
		the same run of the pipeline. If you wish to skip Step 1 of the pipeline and run BWA immediately, all input Fastq files
		must be decompressed, merged (if previously split), and end in .fastq (not .fq or .gz). If starting at a later step in
		the pipeline, SAM files must end in .sam, BAM files must end in .bam, and BAM index files must end in .bai.

		Fastq-SampleID File:
		When working with Fastq files, you must provide a tab-delimited file that lists the SampleIDs corresponding to each Fastq
		file. You only need to list the R1 files, but it's okay to include R2 files, as they'll be ignored. The first column of the 
		SampleID file should hold the filenames of the Fastq files (without any paths), and the second column should contain the
		SampleIDs. Technically, split Fastq files only need to have one of their constituent files included in the list,
		but it's fine to list them all.

		Dealing with Samples with More Than One BAM File:
		When multiple BAM files share a SampleID (this is determined by examining the SM field of the BAM headers), they are merged together using
		GATK's PrintReads tool. Then, the merged file undergoes the Mark Duplicates and Indel Realignment steps again.

	EOF
	$message =~ s/\t//g;
	print $message;
}
