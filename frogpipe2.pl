#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use Pod::Usage;

my $help = 0;
my $QConly = 0;
my $R1File;
my $R2File;
my $adapters;
my $outDir;
my $name = "SampleX";
my $CPU = 6;

GetOptions  ("forward|R1=s"      => \$R1File,
             "reverse|R2=s"      => \$R2File,
             "adapters=s"        => \$adapters,
             "name=s"            => \$name,
             "outdir=s"          => \$outDir,
             "qc_only"           => \$QConly,
             "threads|CPU=i"     => \$CPU,
             "help|man"          => \$help) || pod2usage(2);

if (!$R1File or !$R2File or !$adapters or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

if (-d $outDir) {
    die "Directory $outDir already exists! Stopping script.\n"
}
mkdir $outDir;
$name = $outDir . "/" . $name;



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
system("Trinity.pl --seqType fq --JM 30G --left $TrinityR1 --right $fqj[1] --CPU $CPU");








# SUBROUTINES
# Scythe for adapter contamination
sub Scythe {
    my %args = @_;
    my($filename, $directories, $suffix) = fileparse($args{fastq}, qr/\.fastq/);
    my $out = $outDir . "/" . $filename . "_scythed.fastq";
    system("scythe -q sanger -a $args{adapters} -o $out $args{fastq}");
    return $out;
}

# Sickle for read trimming
sub Sickle {
    my @out;
    push(@out, ($name . "R1_sic"));
    push(@out, ($name . "R2_sic"));
    push(@out, ($name . "sing_sic"));
    system("sickle pe -t sanger -n -f $R1Scythe -r $R2Scythe -o $out[0] -p $out[1] -s $out[2]");
    return @out;
}

# Fastq-join for overlapping reads
sub FastqJoin {
    my $fqj_name = $name . "_scy_sic_fqj_%.fastq";
    my @out;
    push(@out,($name . "_scy_sic_fqj_un1.fastq"));
    push(@out,($name . "_scy_sic_fqj_un2.fastq"));
    push(@out,($name . "_scy_sic_fqj_join.fastq"));
    system("fastq-join -v ' ' $sickle[0] $sickle[1] -o $fqj_name");
    return @out; # this contains the filenames created by fastq-join
}

sub PrepareReadsForTrinity {
    my $join_singles = $name . "_scy_sic_joined_and_singles.fastq";
    my $R1_and_joinsingles = $name . "_scy_sic_R1_and_joinedsingles.fastq";
    system("cat $fqj[2] . $sickle[2] > $join_singles"); #Combine joined and sickle singletons
    system("cat $fqj[0] . $join_singles > $R1_and_joinsingles"); #add all single reads to end of R1 file
    return($R1_and_joinsingles);    
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