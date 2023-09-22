#!/usr/bin/perl -w

=head1 Description

	Merge isolate mutation profile to unique mutation analysis
	
Usage

	perl strainVar.pl -WorkPath <path> -seqid <sequence id> -lineage <lineage information> -varProfile <comparative-results.tsv>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-seqid	[str]	Input the sequence id
	
	-lineage	[str]	Input the lineage information
	
	-varProfile	[str] Input the comparative-results.tsv
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.05.51 15:22 0.0.2
	
=cut

use warnings;
use strict;
use Getopt::Long;

my ($analysis_path, $seqid, $lineage, $varfile, $help) = ("", "", "", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'seqid=s' => \$seqid,
	'lineage=s' => \$lineage,
	'varProfile=s' => \$varfile,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$seqid) or (!$lineage) or (!$varfile) or ($help));

my @out;
push @out, join("\t", ("Seq ID", "Lineage", "Mutation type", "Site in Gene", "Site in Genome", "Mutation", "Gene", "Amino acid change", "Synonymous or non-syn"))."\n";

open(FILE, "<", $analysis_path.$varfile);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		if($line[3] eq "Non-conding region"){
			if($line[1] eq "SNP"){
				push @out, join("\t", ($seqid, $lineage, $line[1], "NA", $line[0], $line[2], $line[3], $line[4], $line[$#line]))."\n";
			}
			elsif($line[1] eq "deletion"){
				my @info = split/\-/,$line[2];
				push @out, join("\t", ($seqid, $lineage, $line[1], "NA", $line[0], "del:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
			}
			elsif($line[1] eq "insertion"){
				my @info = split/\-/,$line[2];
				push @out, join("\t", ($seqid, $lineage, $line[1], "NA", $line[0], "ins:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
			}
		}
		else{
			if($line[2] =~ m/\(/){
				my @res = split/\(/,$line[2];
				$res[1] =~ s/\)//;
				if($line[1] eq "SNP"){
					my @SNP = split/\d+/,$res[1];
					my @pos = split/\D/,$res[1];
					push @out, join("\t", ($seqid, $lineage, $line[1], $pos[1], $line[0], $SNP[0]."->".$SNP[1], $line[3], $line[4], $line[$#line]))."\n";
				}
				elsif($line[1] eq "deletion"){
					my @info = split/\-/,$res[1];
					push @out, join("\t", ($seqid, $lineage, $line[1], $info[0]."~".$info[2], $line[0], "del:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
				}
				elsif($line[1] eq "insertion"){
					my @info = split/\-/,$res[1];
					push @out, join("\t", ($seqid, $lineage, $line[1], $info[0], $line[0], "ins:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
				}
			}
			else{
				if($line[1] eq "SNP"){
					my @SNP = split/\d+/,$line[2];
					my @pos = split/\D/,$line[2];
					push @out, join("\t", ($seqid, $lineage, $line[1], $pos[1], $line[0], $SNP[0]."->".$SNP[1], $line[3], $line[4], $line[$#line]))."\n";
				}
				elsif($line[1] eq "deletion"){
					my @info = split/\-/,$line[2];
					$info[0] =~ s/del\://;
					push @out, join("\t", ($seqid, $lineage, $line[1], $info[0]."~".$info[2], $line[0], "del:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
				}
				elsif($line[1] eq "insertion"){
					my @info = split/\-/,$line[2];
					$info[0] =~ s/ins\://;
					push @out, join("\t", ($seqid, $lineage, $line[1], $info[0], $line[0], "ins:".$info[1], $line[3], $line[4], $line[$#line]))."\n";
				}
			}
		}
	}
}
close FILE;

open(OUT, ">", $analysis_path."Strain.var");
print OUT @out;
close OUT;
