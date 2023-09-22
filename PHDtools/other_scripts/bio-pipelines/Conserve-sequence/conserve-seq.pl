#!/usr/bin/perl -w

=head1 Description
	Find the most conserved sequence of a special pathogen
Usage
	perl conserve-seq.pl -WorkPath <path> -reference <reference genome> -isolates <isolate genomes> -algo <max or shannon> -targetLen <target length> -StorePath <output path>
Parameters
	-WorkPath	[str]	Input the analysis path
	-reference	[str]	Input the reference genome
	-isolates	[str]	Input the isolate genomes
	-algo	[str]	Choose one of the method to calculate the conservatism score [max, shannon]
	-targetLen	[int]	Input the target length
	-StorePath	[str]	Input the output path
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.16 21:11 0.0.3
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $conserveregion = "/var/www/other_scripts/bio-pipelines/Conserve-sequence/Functions/conserveregion.pl";
my $findconserve = "/var/www/other_scripts/bio-pipelines/Conserve-sequence/Functions/findconserve.pl";

my $analysis_path = "";
my ($isolates, $reference, $algorithm, $targetLen, $help) = ("", "", "", 200, "");
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'isolates=s' => \$isolates,
	'reference=s' => \$reference,
	'algo=s' => \$algorithm,
	'targetLen=i' => \$targetLen,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$isolates) or (!$reference) or (!$algorithm) or (!$targetLen) or (!$storage_path) or ($help));

my @shell;
push @shell, "#!/usr/bin/bash\n";
my $refseq = Bio::SeqIO -> new(
	-file => $analysis_path.$reference,
	-format => 'fasta',
);
my $refgenome = $refseq -> next_seq;
my $sequence = $refgenome -> seq;
$sequence = ">Ref\n".$sequence."\n";
unlink $analysis_path.$reference;
open(OUT, ">", $analysis_path.$reference);
print OUT $sequence;
close OUT;
push @shell, $perl." ".$conserveregion." -WorkPath ".$analysis_path." -reference ".$analysis_path.$reference." -isolates ".$analysis_path.$isolates." -algo ".$algorithm."\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $perl." ".$findconserve." -WorkPath ".$analysis_path." -basescore ".$analysis_path."base_score.tsv -stats ".$analysis_path."stat.tsv -targetLen ".$targetLen."\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "zip -qjr ".$analysis_path."MostConserve_info.zip ".$analysis_path."Conserved_region_info.tsv ".$analysis_path."conserved-sequence.fasta\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."MostConserve_info.zip ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
