#!/usr/bin/perl -w

=head1 Description
	Tidy the assembled sequences
Usage
	perl tidyctgs.pl -WorkPath <analysis path> -contigs <Final assembly contigs>
Parameters
	-WorkPath	[str]	Input the analysis path
	-contigs	[str]	Input the final assembly contigs
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.08.25 13:56 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my ($analysis_path, $assemble, $help) = ("", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'contigs=s' => \$assemble,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$assemble) or ($help));

my @out;
my $num = 1;
my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$assemble, -format => 'fasta');
while(my $seq = $fasta -> next_seq){
	my $seq_seq = $seq -> seq;
	push @out, ">Contigs".$num."|".length($seq_seq)."bp|\n".$seq_seq."\n";
	$num += 1;
}
open(OUT, ">", $analysis_path."Final_assembly.fasta");
print OUT @out;
close OUT;
