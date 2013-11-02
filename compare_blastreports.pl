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
             "targets=s"           => \$targets,
             "assembly=s"          => \$assembledContigs,
             "reducedtargets=s"    => \$reducedTargets,
             "reducedassembly=s"   => \$reducedAssembly) || die "couldn't retrieve arguments with GetOpt::Long.\n";

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

my %blastHash;
my %matchesHash;
my $targs_2_assembly = Bio::SearchIO->new(-file => $blast_targets_to_assembly,
                                          -format => 'blast');
while (my $result = $targs_2_assembly->next_result() ) {
    my $name = $result->query_name();
    my @hitArray;
    while (my $hit = $result->next_hit) {
        push(@hitArray,$hit->name());
        $matchesHash{$hit->name()}++;
    }
    $blastHash{$name} = \@hitArray; # number of keys would be scalar(@{$blastHash{$name}})
}


my $targetNew = Bio::SeqIO->new(-file => ">$reducedTargets",
                                -format => "fasta");
my $assemblyNew = Bio::SeqIO->new(-file => ">$reducedAssembly",
                                -format => "fasta");

foreach my $key (sort keys %blastHash) {
    if (scalar(@{$blastHash{$key}}) == 1 and $matchesHash{${$blastHash{$key}}[0]} == 1) {
        $targetNew->write_seq($targetHash{$key});
        #print ${$blastHash{$key}}[0] . "\n";
        $assemblyNew->write_seq($assemblyHash{${$blastHash{$key}}[0]});
    }    
}
