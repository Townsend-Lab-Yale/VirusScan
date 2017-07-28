#!/usr/bin/env perl
use strict;
use warnings;
#use POSIX;

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m"; 

my $version = "sgg0";

#usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the virus discovery pipeline on slurm cluster.
Pipeline version: $version
$yellow		Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green		[1]  Run bwa
$red		[2, or <=22]  Split files for RepeatMasker
			[3 or <=23]  Submit RepeatMasker job array
$yellow		[4 or <=24]  Sequence Qulity Control
$green		[5 or <=25]  Split files for Blast Reference Genome
			[6 or <=26]  Submit Blast Reference Genome job array
			[7 or <=27]  Parse Reference Genome Blast result
$gray		[8 or <=28]  Pool and split files for BlastN
			[9 or <=29]  Submit BlastN job array
			[10 or <=30] Parse BlastN result
			[11 or <=31] Get summary of BlastN
$purple		[12 or <=32] Assignment report for each sample
			[13 or <=33] Assignment summary for each sample
			[14 or <=34] Generate report for the run
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
	$run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 17) || ($step_number >= 22));


#####################################################################################
# values need to be modified to adapt to local environment
my $HOME = $ENV{HOME};
open(my $fwdFile, $HOME."/.forward");
my $email = <$fwdFile>;
my $job_rex = qr/\d+/;  # to match job id from sbatch output (e.g. "Submitted batch job 2653294")

# software path
#my $cd_hit = "/gscuser/mboolcha/software/cdhit/cd-hit-est";
my $repeat_masker = "RepeatMasker";
my $blastn = "blastn";
#my $blastx = "/gscuser/scao/tools/software/ncbi-blast+/bin/blastx";

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# path and name of databases
#my $db_BN = "/gscuser/scao/gc3027/nt/nt";
#my $db_BX = "/gscuser/scao/gc3027/nr/nr";
#my $bwa_ref = "/gscuser/scao/gc3027/fasta/virus/virusdb_082414.fa";

my $db_BN = "nt";
my $db_BX = "nr";
my $genomes_dir = $ENV{'VIRUSSCAN_GENOMES_DIR'};
my $bwa_ref = $genomes_dir."/nt012414_virus_abbr_cdhit98.fa";

# reference genome taxonomy classification and database location.
# It's better to change $refrence_genome_taxonomy and $reference_genome based on the data being analyzed.
my $refrence_genome_taxonomy = "";
my $reference_genome = "";

#if ($ref_genome_choice == 1) {
#	$refrence_genome_taxonomy = "Homo"; # use Bacteria, Homo, Phage, Fungi, Mus, other

	# path to the reference genome	
#	$reference_genome = "/gscmnt/gc3027/dinglab/medseq/human70.37/humandnacdna.fa";
#}

$refrence_genome_taxonomy = "Homo";

$reference_genome = $genomes_dir."/humandnacdna";  # orig "/.../medseq/human70.37/humandnacdna.fa"

#####################################################################################
# everything else below should be automated
my $working_name= (split(/\//,$run_dir))[-2];

# To run jobs faster, split large fasta files to small ones. Split to specific number of 
# files instead of specific sequences in each small file, because the number of job array 
# cannot be determined if spliting to specific number of sequences in each file. Job 
# number is required by qsub ${SGE_TASK_ID}. The minimum size of each file is 4kb. 
# The number of files should be determined accourding to CPUs available in the computer
# cluster.

# The number of small fasta files to split to from a large file for RepeatMasker
my $file_number_of_RepeatMasker = 100; #default 
# the number of small fasta files to split to from a large file for Blast_Reference_Genome
my $file_number_of_Blast_Ref_Genome = 100; #default
# the number of small fasta files to split to from a large file for Blast_N
my $file_number_of_Blast_N = 10; #default. Require simultaneous db access.
# the number of small fasta files to split to from a large file for Blast_X
#my $file_number_of_Blast_X = 200; #default

#store job files here
my $HOME1=$HOME;

#store job files here
if (! -d $HOME1."/tmp") {
    `mkdir $HOME1"/tmp"`;
}
my $job_files_dir = $HOME1."/tmp";

#store SGE output and error files here
if (! -d $HOME1."/SGE_DIR") {
    `mkdir $HOME1"/SGE_DIR"`;
}
my $lsf_file_dir = $HOME1."/SGE_DIR";

# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/env perl ".$run_script_path."/";

my $hold_RM_job = "norm"; 
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

#directory suffix constants
my $REPEAT_MASKER_DIR_SUFFIX = "fa.cdhit_out_RepeatMasker";
my $BLAST_RefG_DIR_SUFFIX = "fa.cdhit_out.masked.goodSeq_RefGblast";
my $BLAST_NT_DIR_SUFFIX = "RefGfiltered_BLASTN";
my $BLASTX_NR_DIR_SUFFIX = "BNFiltered_BLASTX_NR";

# get sample list in the run, name should not contain "."
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
&check_input_dir($run_dir);

# start data processsing
if ($step_number < 14 || $step_number>=22) {
	#begin to process each sample
	for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_dir_list[$i];
		if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
				$current_job_file="";
				if ($step_number == 0 || $step_number>=22) {#run the whole pipeline
					######################################################################
					#cd-hit
                    if($step_number==0) 
					{ &bsub_bwa();}

					######################################################################
					#RepeatMasker
					#split file for RepeatMasker
					#my $f_fa=$sample_full_path.".fa"; 
					#if(! -f $f_fa) { next; }
	
					if($step_number<=22) 
					{
					 &split_for_RepeatMasker(); }
                                        
					#submit RepeatMasker job array
					if($step_number<=23) 
					{
					&submit_job_array_RM();
					$hold_RM_job=$current_job_file; # to limit number repeatmasker jobs run in the cluster at the same time. Can be removed if the cluster is able to handle the volumn of data input/output. 
                                        }
					######################################################################
					#Sequence Quality Control
					if($step_number<=24)
					{ &seq_QC();}

					######################################################################
					#BLASTn against Reference Genome
					if($step_number<=25) 
					{
					&split_for_blast_RefG();}

					#submit Blast RefG job array
					if($step_number<=26) 
					{
					&submit_job_array_blast_RefG();}

					if($step_number<=27)
					{
					#parser Blast RefG file
					&parse_blast_RefG();}


					######################################################################
					#BLASTn against nt
					#pool and split files for BLASTn
					if($step_number<=28)
					{
					&pool_split_for_blast_N();}

					#submit BLASTn job array
					if($step_number<=29)
					{
					&submit_job_array_blast_N();}

					#parser BLASTn output file
					if($step_number<=30)
					{
					 &parse_blast_N();}

                    if($step_number<=31)
                    {
                     &blast_S();}

					if($step_number<=32){
					&report_for_each_sample();}

					#Assignment summary for each sample
					if($step_number<=33) {
					&summary_for_each_sample();}
				
				######################################################################
				#run the pipeline step by step
				}elsif ($step_number == 1) {
					&bsub_bwa();
				}elsif ($step_number == 2) {
					&split_for_RepeatMasker(1);
				}elsif ($step_number == 3) {
					&submit_job_array_RM(1); 
					$hold_RM_job=$current_job_file; # to limit number of repeatmasker jobs
				}elsif ($step_number == 4) {
					&seq_QC(1);
				}elsif ($step_number == 5) {
					&split_for_blast_RefG(1);
				}elsif ($step_number == 6) {
					&submit_job_array_blast_RefG(1);
				}elsif ($step_number == 7) {
					&parse_blast_RefG(1);
				}elsif ($step_number == 8) {
					&pool_split_for_blast_N(1);
				}elsif ($step_number == 9) {
					&submit_job_array_blast_N(1);
				}elsif ($step_number == 10) {
					&parse_blast_N(1);
				}elsif ($step_number == 11) {
                                        &blast_S(1);
				}elsif ($step_number == 12) {
					 &report_for_each_sample(1);
				}elsif ($step_number == 13) {
					 &summary_for_each_sample(1);
				}
			}
		}
	}
}

##########################################################################################
# generate report for the run
if (($step_number == 0) || ($step_number == 14) || ($step_number>=22)) {
	# my $step_by_step = 0;
	# if ( $step_number == 14 ){
	# 	$step_by_step = 1;
	# }
	print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
	$current_job_file = "Run_report_".$$.".sh"; 
	open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
	print REPRUN "#!/bin/bash\n";

    print REPRUN "#SBATCH -c 1\n";
    print REPRUN "#SBATCH --mem 40G\n";
    #print REPRUN "#SBATCH -p ding-lab\n";
	print REPRUN "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print REPRUN "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print REPRUN "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print REPRUN "#SBATCH -d afterok:".$hold_job_file."\n"; }
	
	print REPRUN "BAD_SEQ=fa.cdhit_out.masked.badSeq\n"; #output of RepeatMasker
	print REPRUN "BAD_SEQ=fa.cdhit_out.masked.badSeq\n"; #output of RepeatMasker

	print REPRUN "OUTPUT=".$run_dir."/Analysis_Report_gi_".$working_name."\n";

	print REPRUN 'if [ -f $OUTPUT ] ',"\n"; # file exist 
	print REPRUN "then\n";
	print REPRUN '	grep "# Finished" ${OUTPUT}',"\n";  
	print REPRUN '	CHECK=$?',"\n";
	print REPRUN '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, file not finish
	print REPRUN "	do\n";
	print REPRUN "		".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
	print REPRUN '		grep "# Finished" ${OUTPUT}',"\n";
	print REPRUN '		CHECK=$?',"\n";
	print REPRUN "	done\n";
	print REPRUN "else\n"; # file does not exist
	print REPRUN "	".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
	print REPRUN '	grep "# Finished" ${OUTPUT}',"\n";
	print REPRUN '	CHECK=$?',"\n";
	print REPRUN '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, file not finish
	print REPRUN "	do\n";
	print REPRUN "		".$run_script_path."generate_final_report_gi.pl ".$run_dir." ".$version,"\n";
	print REPRUN '		grep "# Finished" ${OUTPUT}',"\n";
	print REPRUN '		CHECK=$?',"\n";
	print REPRUN "	done\n";
	print REPRUN "fi\n";
	close REPRUN;
	close REPRUN;
    $bsub_com = "sbatch  $job_files_dir/$current_job_file\n";
	#$bsub_com = "qsub -V -P long -hold_jid $working_name -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;

}

#######################################################################
# send email to notify the finish of the analysis
if (($step_number == 0) || ($step_number == 15) || ($step_number>=22)) {
	print $yellow, "Submitting the job for sending an email when the run finishes ",$sample_name, "...",$normal, "\n";
	$current_job_file = "Email_run_".$$.".sh";
	open(EMAIL, ">$job_files_dir/$current_job_file") or die $!;
	print EMAIL "#!/bin/bash\n";
    print EMAIL "#SBATCH -c 1\n";
    print EMAIL "#SBATCH -o $lsf_file_dir","\n";
    print EMAIL "#SBATCH -e $lsf_file_dir","\n";
    print EMAIL "#SBATCH -J $current_job_file\n";
	print EMAIL "#SBATCH -d afterok:".$hold_job_file."\n";
	print EMAIL $run_script_path."send_email.pl ".$run_dir." ".$email."\n";
	close EMAIL;
    $bsub_com = "sbatch  $job_files_dir/$current_job_file\n";
	#$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}
#######################################################################
if ($step_number == 0) {
	print $green, "All jobs are submitted! You will get email notification when this run is completed.\n",$normal;
}

exit;


########################################################################
# subroutines 

sub check_input_dir {
	my ($input_dir) = @_;
	my $have_input_sample = 0;
	
	# get sample list in the run, name should not contain "."
	opendir(DH, $input_dir) or die "Cannot open dir $input_dir: $!\n";
	my @sample_list = readdir DH;
	close DH;
	
	for (my $i=0;$i<@sample_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_list[$i];
		if (!($sample_name =~ /\./)&&!($sample_name =~/Analysis_/)) {
			$have_input_sample = 1;
			$sample_full_path = $input_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				my $input_file = $input_dir."/".$sample_name."/".$sample_name.".bam";
				if (!(-e $input_file)) { # input file does not exist
					print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
					die;
				}
			}
			else { # input sample directory does not exist
				print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
				die;
			}
		}
	}

	if (!($have_input_sample)) { # does not have any input sample directory
		print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
		die;
	}

}

########################################################################
########################################################################
sub bsub_bwa{

    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    $current_job_file = "j1_bwa_".$sample_name.$$.".sh";

    my $IN_bam = $sample_full_path."/".$sample_name.".bam";

    if (! -e $IN_bam) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
    }

    open(BWA, ">$job_files_dir/$current_job_file") or die $!;
    print BWA "#!/bin/bash\n";
    print BWA "#SBATCH -c 1\n";
    print BWA "#SBATCH --mem 20G\n";
    print BWA "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print BWA "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print BWA "#SBATCH -J $current_job_file\n";
    print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
    print BWA "BWA_sai=".$sample_full_path."/".$sample_name.".sai\n";
    #print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".sam\n";
    #print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
    #print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
    print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
    print BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
	#print BWA 
	print BWA 'if [ ! -s $BWA_mapped ]',"\n";
    print BWA "    then\n";
	print BWA "rm \${BWA_sai}","\n";
	print BWA "rm \${BWA_fq}","\n";
	#print BWA "mkfifo \${BWA_sai}","\n";
	print BWA "mkfifo \${BWA_fq}","\n";
	#0x100: secondary alignment
	#0x800: supplementary alignment
    #H: Hard clipping
	#S: Soft clipping
	print BWA "samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}\' | samtools view -u - | bedtools bamtofastq -i - -fq \${BWA_fq} \&","\n";
    #print BWA "bwa aln $bwa_ref -b0 \${BWA_IN} > \${BWA_sai} \&","\n";	
    print BWA "bwa aln $bwa_ref \${BWA_fq} > \${BWA_sai}","\n";
    print BWA 'rm ${BWA_fq}',"\n";
	print BWA "mkfifo \${BWA_fq}","\n";
    print BWA "samtools view -h \${BWA_IN} | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); \$flag=\$ss[1]; \$cigar=\$ss[5]; if(\$ss[0]=~/^\@/ || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && ((\$flag & 0x4) || (\$cigar=~/S/))) || (!((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) && (\$ss[2]=~/^gi/))) { print \$line;}\' | samtools view -u - | bedtools bamtofastq -i - -fq \${BWA_fq} \&","\n";
	#print BWA "samtools view -h \${BWA_IN} | gawk \'{if (substr(\$1,1,1)==\"\@\" || (and(\$2,0x4) || and(\$2,0x8) )) print}\' | samtools view -Sb - | bamtools convert -format fastq > \${BWA_fq} \&","\n";
	print BWA "bwa samse $bwa_ref \${BWA_sai} \${BWA_fq} | grep -v \@SQ | perl -ne \'\$line=\$_; \@ss=split(\"\\t\",\$line); if(\$ss[2]=~/^gi/) { print \$line; }\' > \${BWA_mapped}","\n";
	print BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    print BWA " 	".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
	print BWA 'rm ${BWA_sai}',"\n";
	print BWA 'rm ${BWA_fq}',"\n";
	print BWA "else\n";
    print BWA "     ".$run_script_path."get_fasta_from_bam_filter.pl \${BWA_mapped} \${BWA_fa}\n";
    print BWA "     ".$run_script_path."trim_readid.pl \${BWA_fa} \${BWA_fa}.cdhit_out\n";
    print BWA "   fi\n";
    close BWA;
    $bsub_com = "sbatch  $job_files_dir/$current_job_file\n";
    $_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub split_for_RepeatMasker {
	#split file for RepeatMasker
	# my ($step_by_step) = @_;
	print STDOUT "Building job 2: split_for_RepeatMasker\n";
	$current_job_file = "j2_".$sample_name."_RM_split_".$$.".sh";
	open(RMSPLIT, ">$job_files_dir/$current_job_file") or die $!;
	print RMSPLIT "#!/bin/bash\n";
	print RMSPLIT "#SBATCH -c 1\n";
    #print RMSPLIT "#SBATCH -p ding-lab\n";
    print RMSPLIT "#SBATCH --mem 10G\n";
    print RMSPLIT "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print RMSPLIT "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print RMSPLIT "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print RMSPLIT "#SBATCH -d afterok:".$hold_job_file."\n"; }
	#####################
	print RMSPLIT "RMSPLIT_IN=".$sample_full_path."/".$sample_name.".fa\n";
	print RMSPLIT "RM_DIR=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
	print RMSPLIT "SAMPLE_DIR=".$sample_full_path."\n\n";
	print RMSPLIT "if [ ! -d \${RM_DIR} ]\n";
	print RMSPLIT "then\n";
	print RMSPLIT "	mkdir \${RM_DIR}\n";
	print RMSPLIT "	".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
	print RMSPLIT "	".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
	print RMSPLIT '	CHECK=$?',"\n";
	print RMSPLIT '	while [ ${CHECK} -eq 10 ]',"\n"; # 10 is the error exit code of check_split_cdhit.pl. It will check whether split_cdhit is correctly completed, if not correctly completed
	print RMSPLIT "	do\n"; # run split and check again
	print RMSPLIT "		".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
	print RMSPLIT "		".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
	print RMSPLIT '		CHECK=$?',"\n";
	print RMSPLIT "	done\n";
	print RMSPLIT "else\n"; # RepeatMasker directory already existed (file already splited)
 	print RMSPLIT "	".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
	print RMSPLIT '	CHECK=$?',"\n";
    #check if spliting file is correctly completed, if not correctly completed. check again
 	print RMSPLIT '	while [ ${CHECK} -eq 10 ]',"\n";
	print RMSPLIT "	do\n";# check again
	print RMSPLIT "		".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out -o \${RM_DIR} -n $file_number_of_RepeatMasker -p ".$sample_name.".fa.cdhit_out_file\n";
	print RMSPLIT "		".$run_script_path."check_split_cdhit.pl \${SAMPLE_DIR}\n";
	print RMSPLIT '		CHECK=$?',"\n";
	print RMSPLIT "	done\n";
	print RMSPLIT "fi\n";
	close RMSPLIT;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";	
	#$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub submit_job_array_RM {
	#submit RepeatMasker job array
	# my ($step_by_step) = @_;
	print STDOUT "Building job 3: submit_job_array_RM\n";
	$current_job_file = "j3_".$sample_name."_RM_".$$.".sh";
	open (RM, ">$job_files_dir/$current_job_file") or die $!;
	print RM "#!/bin/bash\n";
	print RM "#SBATCH -c 1\n";
	#print RM "#SBATCH -p ding-lab\n";
    print RM "#SBATCH --mem 10G\n";
    print RM "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print RM "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print RM "#SBATCH -J $current_job_file\n";
	print RM "#SBATCH --array=1-$file_number_of_RepeatMasker\n";
	if ($hold_job_file) { print RM "#SBATCH -d afterok:".$hold_job_file."\n"; }
	print RM "RM_IN=".$sample_full_path."/".$sample_name.".fa\n";
	#####################
	print RM "RM_dir=".$sample_full_path."/".$sample_name.".$REPEAT_MASKER_DIR_SUFFIX\n";
	#print RM "#\$ -t 1-$file_number_of_RepeatMasker:1","\n";
	print RM "RMOUT=",'${RM_dir}',"/".$sample_name.".fa.cdhit_out_file".'${SLURM_ARRAY_TASK_ID}'.".fa.masked","\n";
	print RM "RMIN=",'${RM_dir}',"/".$sample_name.".fa.cdhit_out_file".'${SLURM_ARRAY_TASK_ID}',".fa\n";
	print RM "RMOTHER=",'${RM_dir}',"/".$sample_name.".fa.cdhit_out_file".'${SLURM_ARRAY_TASK_ID}'.".fa.out","\n\n";	
	print RM 'if [ -f $RMIN ]',"\n"; # input file exist
	print RM "then\n";
	print RM '  if [ ! -s $RMOUT ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
    print RM "  then\n";
    #print RM '      while [ ! -s $RMOUT ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
   # print RM "      do\n"; # run RepeatMasker until it finishes
    print RM "          $repeat_masker -pa 4 \$RMIN \n";
   # print RM "      done\n";
    print RM "  fi\n\n";
	print RM '	if [ ! -f $RMOTHER ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
	print RM "	then\n";
	print RM '		while [ ! -f $RMOTHER  ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
	print RM "		do\n"; # run RepeatMasker until it finishes
	print RM "			$repeat_masker -pa 4 \$RMIN \n"; 
	print RM " 		done\n";
	print RM "	fi\n\n";
 	print RM '	if [ ! -f $RMOUT ]',"\n"; #sometimes repeatmasker does not find any repeat in input files, in these cases no .masked file will be generated.
	print RM "	then\n";
	print RM '		cp ${RMIN} ${RMOUT}',"\n";
	print RM "	fi\n";
	print RM "fi\n";
	close RM;
    $bsub_com = "sbatch  $job_files_dir/$current_job_file\n";
	#print $bsub_com, "\n";
        #$bsub_com = "qsub -V -l h_vmem=4G  -hold_jid $hold_job_file,$hold_RM_job -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub seq_QC {
	# my ($step_by_step) = @_;
	print STDOUT "Building job 4: seq_QC\n";
	$current_job_file = "j4_".$sample_name."_QC_".$$.".sh";
	open(QC, ">$job_files_dir/$current_job_file") or die $!;
	print QC "#!/bin/bash\n";
    print QC "#SBATCH -c 1\n";
    print QC "#SBATCH --mem 10G\n";
    print QC "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print QC "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print QC "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print QC "#SBATCH -d afterok:".$hold_job_file."\n"; }
	#####################
	print QC "SAMPLE_DIR=".$sample_full_path."\n";
	print QC "QC_OUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n\n";
	print QC "f_fa=".$sample_full_path."/".$sample_name.".fa\n";
	print QC 'if [ ! -f $QC_OUT ] && [ -s $f_fa ]',"\n";
	print QC "then\n";
	print QC "	".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
	print QC "	".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
	print QC '	CHECK=$?',"\n";
	print QC '	while [ ${CHECK} -eq 10 ]',"\n";#10 is the exit code of check_SequenceQualityControl.pl if it is not correctly completed.
	print QC "	do\n";#run split and check again
	print QC "		".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
	print QC "		".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";	
	print QC '		CHECK=$?',"\n";
 	print QC "	done\n";
	print QC "else\n";
	print QC "	".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
	print QC '	CHECK=$?',"\n";
    #check if parsed file is completed, if not completed. check again
	print QC '	while [ ${CHECK} -eq 10 ]',"\n";
	print QC "	do\n";#run parser again
	print QC "		".$run_script_path."SequenceQualityControl.pl ".$sample_full_path."\n";
	print QC "		".$run_script_path."check_SequenceQualityControl.pl \${SAMPLE_DIR}\n";
	print QC '		CHECK=$?',"\n";
	print QC '              CHECK=1',"\n";
    print QC "	done\n";
	print QC "fi\n";
	close QC;
    $bsub_com = "sbatch  $job_files_dir/$current_job_file\n";
	#$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub split_for_blast_RefG{
	#split file for RefG blast
	# my ($step_by_step) = @_;
	print STDOUT "Building job 5: split_for_blast_RefG\n";
	$current_job_file = "j5_".$sample_name."_RefG_split_".$$.".sh";
	open(RefGS, ">$job_files_dir/$current_job_file") or die $!;
	print RefGS "#!/bin/bash\n";
	print RefGS "#SBATCH -c 1\n";
    print RefGS "#SBATCH --mem 10G\n";
    print RefGS "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print RefGS "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print RefGS "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print RefGS "#SBATCH -d afterok:".$hold_job_file."\n"; }
	############################
	print RefGS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
	print RefGS "SAMPLE_DIR=".$sample_full_path."\n\n";
	print RefGS 'if [ ! -d $RefG_DIR ]',"\n";
	print RefGS "then\n";
	print RefGS "	mkdir \${RefG_DIR}\n";
	print RefGS "	".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
	print RefGS "	".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
	print RefGS '	CHECK=$?',"\n";	
	print RefGS '	while [ ${CHECK} -eq 10 ]',"\n";#10 is the error exit code of it is not correctly completed.
	print RefGS "	do\n";#run split and check again
	print RefGS "		".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
	print RefGS "	".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
	print RefGS '		CHECK=$?',"\n";
	print RefGS "	done\n";
	print RefGS "else\n";
	print RefGS "	".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
	print RefGS '	CHECK=$?',"\n";
    #check if parsed file is completed, if not completed. check again
	print RefGS '	while [ ${CHECK} -eq 10 ]',"\n";
	print RefGS "	do\n";#run parser again
	print RefGS "		".$run_script_path."split_fasta.pl -i ".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq -o \${RefG_DIR} -n $file_number_of_Blast_Ref_Genome -p ".$sample_name.".fa.cdhit_out.masked.goodSeq_file\n";
	print RefGS "		".$run_script_path."check_split_RefG.pl \${SAMPLE_DIR}\n";
	print RefGS '		CHECK=$?',"\n";
	print RefGS "	done\n";
	print RefGS "fi\n";
	close RefGS;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub submit_job_array_blast_RefG{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 6: submit_job_array_blast_RefG\n";
	$current_job_file = "j6_".$sample_name."_BRefG_".$$.".sh";
	open (RefG, ">$job_files_dir/$current_job_file") or die $!;
	print RefG "#!/bin/bash\n";
	print RefG "#SBATCH -c 1\n";
    print RefG "#SBATCH --mem 20G\n";
    print RefG "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print RefG "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print RefG "#SBATCH -J $current_job_file\n";
    print RefG "#SBATCH --array=1-$file_number_of_Blast_Ref_Genome\n";
	if ($hold_job_file) { print RefG "#SBATCH -d afterok:".$hold_job_file."\n"; }
	
	####################
	print RefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
	#print RefG "#\$ -t 1-$file_number_of_Blast_Ref_Genome:1","\n"; #the number must be a digital value in the .sh job file, cannot be calculated when the job submitted
	print RefG "BlastRefGOUT=",'${RefG_DIR}',"/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${SLURM_ARRAY_TASK_ID}',".RefGblast.out\n";
	print RefG "QUERY=",'${RefG_DIR}',"/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${SLURM_ARRAY_TASK_ID}'.".fa\n\n";
	print RefG 'if [ -s $QUERY ]',"\n"; #modified by song: check if a file is empty.
	print RefG "then\n";
	#if blast output file does not exist, do blast and check the completeness of output
	print RefG '	if [ ! -f $BlastRefGOUT ]',"\n";
	print RefG "	then\n";
	print RefG "		$blastn -evalue  1e-9   -show_gis  -num_threads 4 -num_descriptions 2  -num_alignments 2  -query  \${QUERY}  -out \${BlastRefGOUT} -db $reference_genome","\n";
	print RefG '		tail -10 ${BlastRefGOUT}|grep Matrix',"\n";
	print RefG '		CHECK=$?',"\n";
	print RefG '		while [ ${CHECK} -eq 1 ]',"\n";
	print RefG "		do\n";
	print RefG "			$blastn -evalue  1e-9   -show_gis -num_threads 4 -num_descriptions 2  -num_alignments 2  -query  \${QUERY}  -out \${BlastRefGOUT} -db $reference_genome","\n";
	print RefG '			tail -10 ${BlastRefGOUT}|grep Matrix',"\n";
	print RefG '			CHECK=$?',"\n";
	print RefG "		done\n";
	#if blast output file exists, check the completeness of output
	print RefG "	else\n";
	print RefG '		tail -10 ${BlastRefGOUT}|grep Matrix',"\n";
	print RefG '		CHECK=$?',"\n";
	print RefG '		while [ ${CHECK} -eq 1 ]',"\n";
	print RefG "		do\n";
	print RefG "			$blastn -evalue  1e-9   -show_gis -num_threads 4 -num_descriptions 2  -num_alignments 2  -query  \${QUERY}  -out \${BlastRefGOUT} -db $reference_genome","\n";
	print RefG '			tail -10 ${BlastRefGOUT}|grep Matrix',"\n";
	print RefG '			CHECK=$?',"\n";
	print RefG "		done\n";
	print RefG "	fi\n";
	print RefG "fi";
	close RefG;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -l h_vmem=10G -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub parse_blast_RefG{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 7: parse_blast_RefG\n";
   # $current_job_file = "j10_".$sample_name."_PBN_".$$.".sh";
    my $BND=$sample_full_path."/".$sample_name.".".$BLAST_RefG_DIR_SUFFIX;
    #if	
    #my $nn1=`tail $BND/*.out | grep Matrix | wc -l`;
    #my $nn2=`ls $BND/*.out | wc -l`;
    #print $nn1,"\n";
    #print $nn2,"\n";
    #if($nn1 != $nn2) { print "resubmitted blastHG for $sample_name","\n"; &submit_job_array_blast_RefG(1);  }
	#else {
	$current_job_file = "j7_".$sample_name."_PRefG_".$$.".sh";
	open (PRefG, ">$job_files_dir/$current_job_file") or die $!;
	print PRefG "#!/bin/bash\n";
	print PRefG "#SBATCH -c 1\n";
    print PRefG "#SBATCH --mem 10G\n";
    print PRefG "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print PRefG "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print PRefG "#SBATCH -J $current_job_file\n";
    print PRefG "#SBATCH --array=1-$file_number_of_Blast_Ref_Genome\n";
	if ($hold_job_file) { print PRefG "#SBATCH -d afterok:".$hold_job_file."\n"; }
	#################################
	print PRefG "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n";
	#print PRefG "#\$ -t 1-$file_number_of_Blast_Ref_Genome:1","\n";#must be a decimal number
	print PRefG "BlastRefGOUT=${sample_name}.fa.cdhit_out.masked.goodSeq_file".'${SLURM_ARRAY_TASK_ID}',".RefGblast.out\n";#name only, not full path
	print PRefG "BlastRefGIN=",'${RefG_DIR}',"/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${SLURM_ARRAY_TASK_ID}'.".fa\n";#full path
	print PRefG "PARSED=",'${RefG_DIR}',"/".$sample_name.".fa.cdhit_out.masked.goodSeq_file".'${SLURM_ARRAY_TASK_ID}'.".RefGblast.parsed\n\n";
	print PRefG 'if [ -s $BlastRefGIN ]',"\n"; # change -f to -s 
	print PRefG "then\n";
	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print PRefG '	if [ ! -f $PARSED ]',"\n";
	print PRefG "	then\n";
	print PRefG "		".$run_script_path."BLASTn_RefGenome_parser.pl \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy\n";
	#check the completeess of parse
	print PRefG '		tail -5 ${PARSED}|grep Summary',"\n";
	print PRefG '		CHECK=$?',"\n";
	# rerun if not completed
	print PRefG '		while [ ${CHECK} -eq 1 ]',"\n";
	print PRefG "		do\n";#run parse again
	print PRefG "			".$run_script_path."BLASTn_RefGenome_parser.pl  \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy \n";
	#check the completeess of parse
	print PRefG '			tail -5 ${PARSED}|grep Summary',"\n";
	print PRefG '			CHECK=$?',"\n";
	print PRefG "		done\n";
	#if the parsed file exists, check the completeness of the parsed file
	print PRefG "	else\n";
	print PRefG '		tail -5 ${PARSED}|grep Summary',"\n";
	print PRefG '		CHECK=$?',"\n";
	print PRefG '		while [ ${CHECK} -eq 1 ]',"\n"; #not complete
	print PRefG "		do\n";
	print PRefG "			".$run_script_path."BLASTn_RefGenome_parser.pl \${RefG_DIR} \${BlastRefGOUT} $refrence_genome_taxonomy \n";
	print PRefG '			tail -5 ${PARSED}|grep Summary',"\n";
	print PRefG '			CHECK=$?',"\n";
	print PRefG "		done\n";
	print PRefG "	fi\n";
	print PRefG "fi";
	close PRefG;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
	#}
}

#####################################################################################

sub pool_split_for_blast_N{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 8: pool_split_for_blast_N\n";
	$current_job_file = "j8_".$sample_name."_BN_split_".$$.".sh";
	open(BNS, ">$job_files_dir/$current_job_file") or die $!;
	print BNS "#!/bin/bash\n";
	print BNS "#SBATCH -c 1\n";
    print BNS "#SBATCH --mem 10G\n";
    print BNS "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print BNS "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print BNS "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print BNS "#SBATCH -d afterok:".$hold_job_file."\n"; }
	############################
	print BNS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
	print BNS "SAMPLE_DIR=".$sample_full_path."\n";
	print BNS "RefGFiltered_fa=".$sample_full_path."/".$sample_name.".RefGfiltered.fa\n";
	print BNS "RefG_DIR=".$sample_full_path."/".$sample_name.".$BLAST_RefG_DIR_SUFFIX\n\n";
	print BNS 'if [ ! -d $BN_DIR ] ',"\n";
	print BNS "then\n";
	print BNS "	mkdir \${BN_DIR}\n";
	print BNS "fi\n";
	print BNS 'if [ -f $RefGFiltered_fa ] ',"\n";
	print BNS "then\n";
	print BNS "	rm \${RefGFiltered_fa}\n";
	print BNS "fi\n";
	print BNS "cat \${RefG_DIR}/*.RefGfiltered.fa >> \${RefGFiltered_fa}\n";
	print BNS "".$run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
	print BNS 'CHECK=$?',"\n";
	print BNS 'while [ ${CHECK} -eq 10 ]',"\n"; #10 is the exit code of check_split_BN.pl. Check whether it is correctly completed, if not rerun split and check again.
	print BNS "do\n";
	# split to -n number of files, this number should be consistent with 
	# the number of blastn job array submitted bellow
	print BNS "	".$run_script_path."split_fasta.pl -i \${RefGFiltered_fa} -o \${BN_DIR} -n $file_number_of_Blast_N -p ".$sample_name.".RefGfiltered.fa_file\n";
	print BNS "	".$run_script_path."check_split_BN.pl \${SAMPLE_DIR}\n";
	print BNS '	CHECK=$?',"\n";
	print BNS "done\n";
	close BNS;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub submit_job_array_blast_N{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 9: submit_job_array_blast_N\n";
	my $BND=$sample_full_path."/".$sample_name.".".$BLAST_NT_DIR_SUFFIX;
  
    #my $nn1=`tail $BND/*.out | grep Matrix | wc -l`;
    #my $nn2=`ls $BND/*.out | wc -l`;

    #print $nn1,"\n";
    #print $nn2,"\n";

    #if($nn1 != $nn2 || $nn2<200) 
	#{
	$current_job_file = "j9_".$sample_name."_BN_".$$.".sh";
	open (BN, ">$job_files_dir/$current_job_file") or die $!;
	print BN "#!/bin/bash\n";
	print BN "#SBATCH -c 1\n";
    print BN "#SBATCH --mem 40G\n";
    print BN "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print BN "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print BN "#SBATCH -J $current_job_file\[1-$file_number_of_Blast_N\]\n";
    print BN "#SBATCH --array=1-$file_number_of_Blast_N\n";
	if ($hold_job_file) { print BN "#SBATCH -d afterok:".$hold_job_file."\n"; }
	#################################
	print BN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
	#print BN "#\$ -t 1-$file_number_of_Blast_N:1","\n"; #must be a decimal number, the value must be determined when this job file is generated. cannot be a variable
	print BN "BlastNOUT=",'${BN_DIR}',"/",$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.out\n";#full path
	print BN "QUERY=",'${BN_DIR}',"/".$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fa\n\n";
	print BN 'if [ -s $QUERY ]',"\n"; #modified by song. check if the file is empty
	print BN "then\n";
	#if the output file does not exist, run and check the completeness of the output file
	print BN '	if [ ! -f $BlastNOUT ]',"\n";
	print BN "	then\n";
	print BN "		$blastn -evalue  1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN","\n";
	print BN '		tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print BN '		CHECK1=$?',"\n";
	print BN '		grep "no longer exists in database" ${BlastNOUT}',"\n"; # one possible blast error message ( see the end of this script).
	print BN '		CHECK2=$?',"\n";
	print BN '		while [ ${CHECK1} -eq 1 ] || [ ${CHECK2} -eq 0 ]',"\n";
	print BN "		do\n";
	print BN "			$blastn -evalue  1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN","\n";
	print BN '			tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print BN '			CHECK1=$?',"\n";
	print BN '			grep "no longer exists in database" ${BlastNOUT}',"\n";#see the end of this script
	print BN '			CHECK2=$?',"\n";
	print BN "		done\n";
	#if the output file exists, check the completeness of the output file
	print BN "	else\n";
	print BN '		tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print BN '		CHECK1=$?',"\n";
	print BN '		grep "no longer exists in database" ${BlastNOUT}',"\n";# one possible blast error (see the end of this script). 
	print BN '		CHECK2=$?',"\n";
	print BN '		while [ ${CHECK1} -eq 1 ] || [ ${CHECK2} -eq 0 ]',"\n";
	print BN "		do\n";
	print BN "			$blastn -evalue  1e-9 -show_gis -num_threads 4 -query \${QUERY} -out \${BlastNOUT} -db $db_BN","\n";
	print BN '			tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print BN '			CHECK1=$?',"\n";
	print BN '			grep "no longer exists in database" ${BlastNOUT}',"\n";#see the end of this script
	print BN '			CHECK2=$?',"\n";
	print BN "		done\n";
	print BN "	fi\n";
	print BN "fi";
	close BN;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -l h_vmem=10G -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
    #}
}

#####################################################################################

sub parse_blast_N{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 10: parse_blast_N\n";
	$current_job_file = "j10_".$sample_name."_PBN_".$$.".sh";
	#my $BND=$sample_full_path."/".$sample_name.".".$BLAST_NT_DIR_SUFFIX;
	#my $nn1=`tail $BND/*.out | grep Matrix | wc -l`;  
    #my $nn2=`ls $BND/*.out | wc -l`; 
	#print $nn1,"\n";
	#print $nn2,"\n";
	#if($nn1 != $nn2) { print "resubmited blastN for $sample_name","\n"; &submit_job_array_blast_N(1);  }
	#else {
    #exit(2); 	
	open (PBN, ">$job_files_dir/$current_job_file") or die $!;
	print PBN "#!/bin/bash\n";
	print PBN "#SBATCH -c 1\n";
    print PBN "#SBATCH --mem 10G\n";
    print PBN "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print PBN "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print PBN "#SBATCH -J $current_job_file\[1-$file_number_of_Blast_N\]\n";
    print PBN "#SBATCH --array=1-$file_number_of_Blast_N\n";
	if ($hold_job_file) { print PBN "#SBATCH -d afterok:".$hold_job_file."\n"; }
	#################################
	print PBN "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
	#print PBN "#\$ -t 1-$file_number_of_Blast_N:1","\n"; #must be a decimal number when the job file is created, cannot be a variable
	print PBN "BlastNOUT=",$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.out\n";#name only, not full path
	print PBN "BlastNIN=",'${BN_DIR}',"/",$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fa\n";#full path
	print PBN "PARSED=",'${BN_DIR}',"/".$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.parsed\n\n";
	print PBN 'if [ -s $BlastNIN ]',"\n"; #song changed -f to -s; 
	print PBN "then\n";
	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print PBN '	if [ ! -f $PARSED ]',"\n";
	print PBN "	then\n";
	print PBN "		".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
	print PBN "		".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print PBN '		CHECK=$?',"\n";
	#check if parsed file is completed, if not completed. run and check again
	print PBN '		while [ ${CHECK} -eq 10 ]',"\n"; #10 is the error exit code of check_Blast_parsed_file.pl if it's not correctly completed.  
	print PBN "		do\n"; #run parser again
	print PBN "			".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
	print PBN "			".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print PBN '			CHECK=$?',"\n";
	print PBN "		done\n";
	#if the parsed file exists, check the completeness of the parsed file
	print PBN "	else\n";
   #     print PBN "             ".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
	print PBN "		".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print PBN '		CHECK=$?',"\n";
	#check if parsed file is completed. If not correctly completed run and check again
	print PBN '		while [ ${CHECK} -eq 10 ]',"\n";
	print PBN "		do\n";  #run parser again
	print PBN "			".$run_script_path."BLASTn_NT_parser.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNOUT}\n";
	print PBN "			".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print PBN '			CHECK=$?',"\n";
	print PBN "		done\n";
	print PBN "	fi\n";
	print PBN "fi";
	close PBN;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub blast_S{

        # my ($step_by_step) = @_;
        print STDOUT "Building job 11: blast_S\n";
        $current_job_file = "j11_".$sample_name."_blastS_".$$.".sh";
        open (PS, ">$job_files_dir/$current_job_file") or die $!;
        print PS "#!/bin/bash\n";
        print PS "#SBATCH -c 1\n";
        print PS "#SBATCH --mem 10G\n";
        print PS "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
        print PS "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
        print PS "#SBATCH -J $current_job_file\[1-$file_number_of_Blast_N\]\n";
        print PS "#SBATCH --array=1-$file_number_of_Blast_N\n";
        if ($hold_job_file) { print PS "#SBATCH -d afterok:".$hold_job_file."\n"; }
        #################################
        print PS "BN_DIR=".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX\n";
        #print PBN "#\$ -t 1-$file_number_of_Blast_N:1","\n"; #must be a decimal number when the job file is created, cannot be a variable
        print PS "BlastNparsed=",$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.parsed\n";#name only, not full path
        print PS "BlastNIN=",'${BN_DIR}',"/",$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fa\n";#full path
        print PS "OUTPUT=",'${BN_DIR}',"/".$sample_name.".RefGfiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.summary\n\n";
        print PS 'if [ -s $BlastNIN ]',"\n"; #song changed -f to -s; 
        print PS "then\n";
        #if the parsed file does not exist, run parser and check the completeness of the parsed file
        print PS '     if [ ! -f $OUTPUT ]',"\n";
        print PS "     then\n";
        print PS "             ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
 		print PS '             grep "Finished summary" ${OUTPUT}',"\n";
        print PS '             CHECK=$?',"\n";
        #check if parsed file is completed, if not completed. run and check again
        print PS '             while [ ${CHECK} -eq 1 ]',"\n"; #10 is the error exit code of check_Blast_parsed_file.pl if it's not correctly completed.  
        print PS "             do\n"; #run parser again
        print PS "                     ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
        print PS '             	       grep "Finished summary" ${OUTPUT}',"\n";
        print PS '                     CHECK=$?',"\n";
        print PS "             done\n";
        #if the parsed file exists, check the completeness of the parsed file
        print PS "     else\n";
        #print PS "             ".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
 		print PS '             grep "Finished summary" ${OUTPUT}',"\n";
        print PS '             CHECK=$?',"\n";
        #check if parsed file is completed. If not correctly completed run and check again
        print PS '             while [ ${CHECK} -eq 1 ]',"\n";
        print PS "             do\n";  #run parser again
        print PS "                     ".$run_script_path."blast_summary.pl ".$sample_full_path."/".$sample_name.".$BLAST_NT_DIR_SUFFIX \${BlastNparsed}\n";
		print PS '                     grep "Finished summary" ${OUTPUT}',"\n";
        print PS '                     CHECK=$?',"\n";
        print PS "             done\n";
        print PS "     fi\n";
        print PS "fi";
        close PS;
        $bsub_com = "sbatch  $job_files_dir/$current_job_file";
        #$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
        $_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
		$hold_job_file = $& if m/$job_rex/;
}


#####################################################################################

sub report_for_each_sample{
	# my ($step_by_step) = @_;
	print STDOUT "Building job 12: report_for_each_sample\n";
	$current_job_file = "j12_".$sample_name."_Rep_".$$.".sh";
	open(REP, ">$job_files_dir/$current_job_file") or die $!;
	print REP "#!/bin/bash\n";
	print REP "#SBATCH -c 1\n";
    print REP "#SBATCH --mem 40G\n";
    print REP "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print REP "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print REP "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print REP "#SBATCH -d afterok:".$hold_job_file."\n"; }
	############################
	print REP "INPUT=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.goodSeq\n";#RepeatMasker QC output
	print REP "REPORT=".$sample_full_path."/".$sample_name.".gi.AssignmentReport\n";
	print REP 'if [ -f $REPORT ] ',"\n"; # report file exist 
	print REP "then\n";
	print REP '	grep "# Finished Assignment Report" ${REPORT}',"\n";  
	print REP '	CHECK=$?',"\n";
	print REP '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, report not finish
	print REP "	do\n";
	print REP "		".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy \n";
	print REP '		grep "# Finished Assignment Report" ${REPORT}',"\n";
	print REP '		CHECK=$?',"\n";
	print REP "	done\n";
	print REP "else\n"; # report file does not exist
	print REP "	".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy \n";
	print REP '	grep "# Finished Assignment Report" ${REPORT}',"\n";  
	print REP '	CHECK=$?',"\n";
	print REP '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, report not finish
	print REP "	do\n";
	print REP "		".$run_script_path."assignment_report_virus_gi.pl ".$sample_full_path." \${INPUT} $refrence_genome_taxonomy \n";
	print REP '		grep "# Finished Assignment Report" ${REPORT}',"\n";
	print REP '		CHECK=$?',"\n";
	print REP "	done\n";
	print REP "fi\n";
	close REP;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -P long -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

#####################################################################################

sub summary_for_each_sample{

	# my ($step_by_step) = @_;
	print STDOUT "Building job 13: summary_for_each_sample\n";
	$current_job_file = "j13_".$sample_name."_Sum_".$$.".sh";

	open(SUM, ">$job_files_dir/$current_job_file") or die $!;
	print SUM "#!/bin/bash\n";
	print SUM "#SBATCH -c 1\n";
    print SUM "#SBATCH --mem 40G\n";
    print SUM "#SBATCH -o $lsf_file_dir","/","$current_job_file.out\n";
    print SUM "#SBATCH -e $lsf_file_dir","/","$current_job_file.err\n";
    print SUM "#SBATCH -J $current_job_file\n";
	if ($hold_job_file) { print SUM "#SBATCH -d afterok:".$hold_job_file."\n"; }
	############################
	print SUM "OUTPUT=".$sample_full_path."/".$sample_name.".gi.AssignmentSummary\n";
	print SUM "BAD_SEQ=".$sample_full_path."/".$sample_name.".fa.cdhit_out.masked.badSeq\n\n"; #output of RepeatMasker
	print SUM 'if [ -f $OUTPUT ] ',"\n"; # summary file exist 
	print SUM "then\n";
	print SUM '	grep "# Finished Assignment Summary" ${OUTPUT}',"\n";  
	print SUM '	CHECK=$?',"\n";
	print SUM '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, file not finish
	print SUM "	do\n";
	print SUM "		".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
	print SUM '		grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
	print SUM '		CHECK=$?',"\n";
	print SUM "	done\n";
	print SUM "else\n"; # file does not exist
	print SUM "	".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
	print SUM '	grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
	print SUM '	CHECK=$?',"\n";
	print SUM '	while [ ${CHECK} -eq 1 ] ',"\n"; # grep unsuccessful, file not finish
	print SUM "	do\n";
	print SUM "		".$run_script_path."assignment_summary_gi.pl ".$sample_full_path." \${BAD_SEQ}\n";
	print SUM '		grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
	print SUM '		CHECK=$?',"\n";
	print SUM "	done\n";
	print SUM "fi\n";
	close SUM;
	$bsub_com = "sbatch  $job_files_dir/$current_job_file";
	#$bsub_com = "qsub -V -P long -N $working_name -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	$_ = qx($bsub_com);  # outputs string of form "Submitted batch job XXXXXXX"
	$hold_job_file = $& if m/$job_rex/;
}

=add
possible blast error
Sequence with id 224967180 no longer exists in database...alignment skipped
Sequence with id 224967180 no longer exists in database...alignment skipped
=cut
