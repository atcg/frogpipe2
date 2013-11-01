#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $R1File;
my $R2File;

GetOptions  ("forward=s"      => \$R1File,
             "reverse=s"      => \$R2File,
             "help|man" => \$help) || pod2usage(2);

if (!$inFile or !$outFile or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}


















#Documentation
__END__

=head1 NAME

perl_template.pl ##CHANGE

=head1 SYNOPSIS 

perl velvetTarget.pl --in <file> --out <file>

 Options:
   -in=s            Input filename
   -out=s           Output filename
   -help|man        Prints out documentation


=head1 DESCRIPTION

Short description of the program.

Can be multiple paragraphs.

=cut