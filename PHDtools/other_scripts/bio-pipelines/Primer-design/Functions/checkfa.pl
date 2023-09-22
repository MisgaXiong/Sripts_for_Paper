#!/usr/bin/perl -w

=head1 Description
	Check the fasta format
Usage
	perl uniq-mutation.pl -WorkPath <path> -target <target sequence file>
Parameters
	-WorkPath	[str]	Input the analysis path
	-target	[str]	Input the target file
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.08.17 20:22 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my ($analysis_path, $target_file, $help) = ("", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'target=s' => \$target_file,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$target_file) or ($help));

my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$target_file, -format => 'fasta');
my $seq = $fasta -> next_seq;
my ($id, $sequence) = ($seq -> id, $seq -> seq);
if(!$id or !$sequence){
	print "FALSE"."\n";
}
else{
	print "TRUE"."\n";
}
