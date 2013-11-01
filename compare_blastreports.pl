#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;

my $blast_targets_to_assembly;
my $blast_assembly_to_targets;
my $targets;
my $assembledContigs;
my $reducedTargets;
my $reducedAssembly;
my $logFile;

GetOptions  ("targets2assembly=s"  => \$blast_targets_to_assembly,
             "assembly2targets=s"  => \$blast_assembly_to_targets,
             "targets=s"           => \$targets,
             "assembly=s"          => \$assembledContigs,
             "reducedtargets=s"    => \$reducedTargets,
             "reducedassembly=s"   => \$reducedAssembly,
             "log=s"               => \$logFile) || die "couldn't retrieve arguments with GetOpt::Long.\n";

# if (!$blast_targets_to_assembly or !$blast_assembly_to_targets or !$targets or !$assembledContigs or !$reducedTargets or !$reducedAssembly or !$logFile) {
#     die "must supply all command line arguments!\n";
# }
# open(my $LOG, ">", $logFile) or die "Couldn't open $logFile for writing.\n";

# Create hashes of all the sequences in the targets and the assembly
my %targetHash;
my $targetIn = Bio::SeqIO->new(-file => "$targets",
                               -format => "fasta");
while (my $seq = $targetIn->next_seq()) { # Create hash of all targets
    # print $seq->display_id() . "\n";
    $targetHash{$seq->display_id()} = $seq;   
}
print "\nFinished putting targets into hash.\n";
my %assemblyHash;
my $assemblyIn = Bio::SeqIO->new(-file=> "$assembledContigs",
                                 -format => "fasta");
while (my $seq = $assemblyIn->next_seq()) { # Create hash of all contigs
    # print $seq->display_id() . "\n";
    $assemblyHash{$seq->display_id()} = $seq;    
}
print "Finished putting assembly contigs into hash.\n";
print scalar(keys %targetHash) . " keys initially in targetHash.\n";
print scalar(keys %assemblyHash) . " keys initially in targetHash.\n";


my $assem2targBLAST2 = Bio::SearchIO->new(-file => "$blast_assembly_to_targets",
                                         -format => "blast");

my %assemblyMatchesHash; # Keys are the names of each contig, values are arrays of the targets they blast match
my %targetsTotalMatchesHash; # Keys are target names, values are number of times the target has been matched by contigs
while (my $result = $assem2targBLAST2->next_result()) {
    my @matches; #this will be an array of target names that are blast hits to each contig
    while (my $hit = $result->next_hit()) {
        push(@matches,$hit->name());
        $targetsTotalMatchesHash{$hit->name()}++;
    }
    $assemblyMatchesHash{$result->query_name()} = \@matches;
}

foreach my $key (sort keys %assemblyHash) {
    #key example: comp9536_c2_seq3
    print "Number of values in matches array: " . scalar(@{$assemblyMatchesHash{$key}}) . "\n";
    if (scalar(@{$assemblyMatchesHash{$key}}) != 1) {
        delete $assemblyHash{$key};    #Get rid of contigs that match multiple targets
    }
}
print scalar(keys %assemblyHash) . " total contigs remaining in assembly after pulling out contigs that don't match a single target.\n";


print scalar(keys %targetsTotalMatchesHash) . " total targets were matched by blasting the assembly to the targets.\n";
my $singlematches=0;
my %singleMatches;
my %multiMatches; # This is a hash where the values are names of targets that were matches for multiple contigs. Values are just "1".
foreach my $key (sort keys %targetsTotalMatchesHash) {
    if ($targetsTotalMatchesHash{$key} == 1) {
        $singlematches++;
        $singleMatches{$key} = 1;
    } else {
        $multiMatches{$key} = 1;
    }
}

print "There are $singlematches total targets that were matched ONCE by a contig.\n";

my $assem2targBLAST3 = Bio::SearchIO->new(-file => "$blast_assembly_to_targets",
                                         -format => "blast");
while (my $result = $assem2targBLAST3->next_result()) {
    while (my $hit = $result->next_hit()) {
        if (exists $multiMatches{$hit->name()}) {
            delete $assemblyHash{$result->query_name()};
        }
    }
}
print scalar(keys %assemblyHash) . " keys in assembly hash after pulling contigs that shared that share the same target match.\n";
print scalar(keys %singleMatches) . " keys in singleMatches hash.\n";


# The name of the matching targets should be in $assemblyMatchesHash{$key}[0]

my $targetNew = Bio::SeqIO->new(-file => ">$reducedTargets",
                                -format => "fasta");
foreach my $key (sort keys %assemblyHash) {
    $targetNew->write_seq($targetHash{$assemblyMatchesHash{$key}[0]});
}

my $assemblyNew = Bio::SeqIO->new(-file => ">$reducedAssembly",
                                -format => "fasta");
foreach my $key (sort keys %assemblyHash) {
    $assemblyNew->write_seq($assemblyHash{$key});
}
