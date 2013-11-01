#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

my $R1File;
my $R2File;
my $adapters;
my $targets;
my $outDir = "frogpipe_run";
my $name = "SampleX";
my $logFile;
my $CPU = 6;
my $QConly = 0;
my $throughAssemblyOnly = 0;
my $help = 0;


GetOptions  ("forward|R1=s"      => \$R1File,
             "reverse|R2=s"      => \$R2File,
             "adapters=s"        => \$adapters,
             "name=s"            => \$name,
             "outdir=s"          => \$outDir,
             "targets=s"         => \$targets,
             "log=s"             => \$logFile,
             "qc_only"           => \$QConly,
             "through_assembly"  => \$throughAssemblyOnly,
             "threads|CPU=i"     => \$CPU,
             "help|man"          => \$help) || pod2usage(2);

if (!$R1File or !$R2File or !$adapters or !$logFile or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

if (-d $outDir) {
    die "Directory $outDir already exists! Stopping script.\n"
}
mkdir $outDir;

open(my $logFH, ">", ($outDir . "/" . $logFile)) or die "Couldn't open log file for writing.\n";

### *** QC STEPS *** ###
#Scythe
my $R1Scythe = Scythe(fastq=>$R1File, adapters=>$adapters);
my $R2Scythe = Scythe(fastq=>$R2File, adapters=>$adapters);

#Sickle
my @sickle = Sickle(); # returns an array of the three sickle output file names in this order: R1, R2, singles
#Fastq-join
my @fqj = FastqJoin(); # returns an array of the three fastq-join output file names in this order: un1, un2, join

# Stop here if you only want QC
if ($QConly) {
    die "Finished with QC steps. Exiting because --qc_only was invoked.\n";
}

### *** ASSEMBLY *** ###
my $TrinityR1 = PrepareReadsForTrinity(); # returns the filename of the combined R1/joined/singles file
my $assembly = Trinity();
# Stop here if you only want to go through assembly
if ($throughAssemblyOnly) {
    die "Finished with Trinity assembly. Exiting because --through_assembly was invoked.\n";
}






# SUBROUTINES
# Scythe for adapter contamination
sub Scythe {
    my $scytheDir = $outDir . "/scythe/";
    unless (-d $scytheDir) {
        mkdir $scytheDir;
    }
    my %args = @_;
    my($filename, $directories, $suffix) = fileparse($args{fastq}, qr/\.fastq/);
    my $out = $scytheDir . $filename . "_scythed.fastq";
    my $scytheOutput = $scytheDir . "Scythelog_" . $filename . ".txt";
    print $logFH "Running scythe on $args{fastq} to find adapter contamination.\n";
    system("scythe -q sanger -a $args{adapters} -o $out $args{fastq} > $scytheOutput 2>&1");
    print $logFH "Finished running scythe on $args{fastq}\n\n";
    return $out;
}

# Sickle for read trimming
sub Sickle {
    my $sickleDir = $outDir . "/sickle/";
    unless (-d $sickleDir) {
        mkdir $sickleDir;
    }
    my @out;
    push(@out, ($sickleDir . $name . "_R1_sic.fastq"));
    push(@out, ($sickleDir . $name . "_R2_sic.fastq"));
    push(@out, ($sickleDir . $name . "_sing_sic.fastq"));
    my $sickleOutput = $sickleDir . "sickleLog.txt";
    print $logFH "Running sickle pe.\n";
    print $logFH "Sickle command: sickle pe -t sanger -n -f $R1Scythe -r $R2Scythe -o $out[0] -p $out[1] -s $out[2]\n";
    system("sickle pe -t sanger -n -f $R1Scythe -r $R2Scythe -o $out[0] -p $out[1] -s $out[2] > $sickleOutput");
    print $logFH "Finished running sickle\n\n";
    return @out;
}

# Fastq-join for overlapping reads
sub FastqJoin {
    my $fqjDir = $outDir . "/fastq-join/";
    unless (-d $fqjDir) {
        mkdir $fqjDir;
    }
    my $fqj_name = $fqjDir . $name . "_scy_sic_fqj_%.fastq";
    my @out;
    push(@out,($fqjDir . $name . "_scy_sic_fqj_un1.fastq"));
    push(@out,($fqjDir . $name . "_scy_sic_fqj_un2.fastq"));
    push(@out,($fqjDir . $name . "_scy_sic_fqj_join.fastq"));
    my $fqjOutput = $fqjDir . "fastq-join.log";
    print $logFH "Running fastq-join on $sickle[0] and $sickle[1]\n";
    system("fastq-join -v ' ' $sickle[0] $sickle[1] -o $fqj_name > $fqjOutput");
    return @out; # this contains the filenames created by fastq-join
}

sub PrepareReadsForTrinity {
    my $assemblyDir = $outDir . "/readsforassembly/";
    unless (-d $assemblyDir) {
        mkdir $assemblyDir;
    }
    my $join_singles = $assemblyDir . $name . "_scy_sic_joined_and_singles.fastq";
    my $R1_and_joinsingles = $assemblyDir . $name . "_scy_sic_R1_and_joinedsingles.fastq";
    system("cat $fqj[2] $sickle[2] > $join_singles"); #Combine joined and sickle singletons
    system("cat $fqj[0] $join_singles > $R1_and_joinsingles"); #add all single reads to end of R1 file
    return($R1_and_joinsingles);    
}

sub Trinity {
    my $trinityDir = $outDir . "/trinity/";
    unless (-d $trinityDir) {
        mkdir $trinityDir;
    }
    my $trinityOutput = $trinityDir . "trinity.log";
    print $logFH "Running Trinity to assembly reads.\n";
    print $logFH "Command: Trinity.pl --seqType fq --JM 16G --output $trinityDir --left $TrinityR1 --right $fqj[1] --CPU $CPU > $trinityOutput\n";
    system("Trinity.pl --seqType fq --JM 16G --output $trinityDir --left $TrinityR1 --right $fqj[1] --CPU $CPU > $trinityOutput");
    print $logFH "Finished running Trinity.\n\n";
    my $assemblyName = $trinityDir . "Trinity.fasta";
    return $assemblyName;
}






#Documentation
__END__

=head1 NAME

frogpipe2.pl

=head1 SYNOPSIS 

perl frogpipe2.pl --R1 <file> --R2 <file> --adapters <file> --targets <file> --name <string> 

 Options:
   --forward|R1=s           R1 reads fastq file
   --reverse|R2=s           R2 reads fastq file
   --adapters=s             Adapter fasta file for trimming purposes
   --name=s                 Sample name to be used as file prefixes
   --targets=s              Fasta file containing targets regions
   --outdir=s               Output directory
   --qc_only                Flag. If included only does scythe/sickle/fastq-join
   --threads                Number of threads used in Trinity
   --help|man               Prints out documentation


=head1 DESCRIPTION

This script runs QC on fastq reads then pipes the reads through the Trinity
assembler.

It then reduces the redundancy of the assembly and does some comparison
between the assembly and the target regions you tried to enrich in your
experiment.

Finally, it maps reads back to the processed assembly and calls SNPs.

=cut
