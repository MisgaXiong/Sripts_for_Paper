#!/usr/bin/perl -w

=head1 Description
	delete fasta sequence from the gff file
=head1 Usage
	perl tidy-gff.pl -WorkPath <path> -gff <gff file>
=head1 Parameters
	-WorkPath		[str]	Analysis path, required
	-gff		[str]	gff file output from prokka, required
	-h/-help		[str]	print help
=head1 Author
	Dongyan Xiong
=head1 Edit Time
	2022.08.20 22:33 0.0.1
=cut

use strict;
use warnings;
use Getopt::Long;

my ($analysis_path, $gff, $help) = ("", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'gff=s' => \$gff,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$gff) or ($help));

my @out;
open(FILE, "<", $analysis_path.$gff);
while(<FILE>){
	chomp;
	if($_ eq "##FASTA"){
        last;
	}
	else{
		push @out, $_."\n";
	}
}
close FILE;
open(OUT, ">", $analysis_path."genome-annotation.gff3");
print OUT @out;
close @out;
