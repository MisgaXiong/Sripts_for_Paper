#!/usr/bin/perl -w

=head1 Description
	Primer design
Usage
	perl primer-design.pl -WorkPath <path> -target <target sequence file> -type <nPCR or qPCR>
	or
	perl primer-design.pl -WorkPath <path> -target <target sequence file> -type <nPCR or qPCR> -consFILE <conserved score file>
Parameters
	-WorkPath	[str]	Input the analysis path
	-target	[str]	Input the target nucleotide sequence but no more than 300 bp.
	-type	[str]	Input the PCR type (nPCR or qPCR)
	-consFILE	[str] Input the converse score file, output from conserved sequence determine
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.08.22 00:53 0.0.1
	2022.07.21 00:57 0.0.2
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $primerdesign = "/var/www/other_scripts/bio-pipelines/Primer-design/Functions/PrimerDesign.pl";

my $analysis_path = "";
my ($target_file, $PCR_type, $cons_file, $help) = ("", "", "", "");
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'target=s' => \$target_file,
	'type=s' => \$PCR_type,
	'consFILE:s' => \$cons_file,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$PCR_type) or (!$target_file) or (!$storage_path) or ($help));


my @shell;
push @shell, "#!/usr/bin/bash\n";
my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$target_file, -format => 'fasta');
my $seq = $fasta -> next_seq;
my $sequence = $seq -> seq;
$sequence = uc($sequence);
die `pod2text $0` if(length($sequence) > 360);
my $checkseq = $sequence;
$checkseq =~ s/A//g;
$checkseq =~ s/T//g;
$checkseq =~ s/C//g;
$checkseq =~ s/G//g;
die `pod2text $0` if(length($checkseq) >= 1);
die `pod2text $0` if($PCR_type ne "qPCR" and $PCR_type ne "nPCR");

if(!$cons_file){
	push @shell, $perl." ".$primerdesign." -WorkPath ".$analysis_path." -sequence ".$sequence." -type ".$PCR_type."\n";
}
else{
	push @shell, $perl." ".$primerdesign." -WorkPath ".$analysis_path." -sequence ".$sequence." -type ".$PCR_type." -consFILE ".$cons_file."\n";
}
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."Primers.tsv ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
