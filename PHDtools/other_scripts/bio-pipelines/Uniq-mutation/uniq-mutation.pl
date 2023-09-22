#!/usr/bin/perl -w

=head1 Description

	Construct the lineage mutation network to identify the unique mutation of some certain lineages
	
Usage

	perl uniq-mutation.pl -WorkPath <path> -reference <reference genome> -isolate <isolate genomes> -lineage <lineage information table> -gff <gff file> -codon <codon table> -p <mutation rate threshold in a lineage> -minNum <minum number of isolate of a lineage> -StorePath <output path>

Parameters

	-WorkPath	[str]	Input the analysis path
	
	-reference	[str]	Input the reference genome
	
	-isolate	[str]	Input the isolate genomes
	
	-lineage	[str]	Input the lineage information table ("Seq ID"	"lineage")
	
	-gff	[str]	Input the gene annotation file
	
	-codon	[int]	Input the translation table
	
	-p	[float]	Input the mutation rate as a threshold in a certain lineage
	
	-minNum	[int]	Input the minum number of isolate in a certain lineage
	
	-StorePath	[str] Input the output path
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.08.04 14:02 0.0.2
	
=cut

use strict;
use warnings;
use Getopt::Long;

my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $batchlin = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/batchlin.pl";
my $lineageVar = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/lineageVarWhole.pl";
my $network = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/lineageWholeVarstat.pl";

my $analysis_path = "";
my ($reference, $isolate, $gff, $lintbl, $help) = ("", "", "", "", "");
my $codon_tbl = 1;
my ($percent, $minNum) = (0.5, 5);
my $storage_path = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'isolate=s' => \$isolate,
	'gff=s' => \$gff,
	'lineage=s' => \$lintbl,
	'codon:i' => \$codon_tbl,
	'p:f' => \$percent,
	'minNum:i' => \$minNum,
	'StorePath=s' => \$storage_path,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$isolate) or (!$gff) or (!$lintbl) or (!$storage_path) or ($help));

my @shell;
push @shell, "#!/usr/bin/bash\n";
push @shell, $perl." ".$batchlin." -WorkPath ".$analysis_path." -reference ".$reference." -isolate ".$isolate." -gff ".$gff." -codon ".$codon_tbl." -lineage ".$lintbl."\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $perl." ".$lineageVar." -WorkPath ".$analysis_path." -TotalVar Total_strainVar.tsv -p ".$percent." -minNum ".$minNum."\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, $perl." ".$network." ".$analysis_path."sig_LineageVar.tsv > ".$analysis_path."lineage-mutation-network.tsv\n";
push @shell, "chmod -R 777 ".$analysis_path."\n";
push @shell, "mkdir ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "mv ".$analysis_path."lineage-mutation-network.tsv ".$storage_path."\n";
push @shell, "chmod -R 777 ".$storage_path."\n";
push @shell, "rm -rf ".$analysis_path."\n";
open(OUT, ">", $analysis_path.".../shell.sh");
print OUT @shell;
close OUT;
