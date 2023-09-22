#!/usr/bin/perl -w

=head1 Description

	Merge the comparative genomics results
	
Usage

	perl mergemutation.pl -WorkPath <path> -snp <snp file> -del <deletion file> -ins <insertion file>
	
	or
	
	perl mergemutation.pl -WorkPath <path> -snp <snp file> -del <deletion file> -ins <insertion file> -AAchange <SNP that change amino acid file> -gff <gff file>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-snp	[str]	Input the SNP file

	-del	[str]	Input the deletion file
	
	-ins	[str]	Input the insertion file
	
	-gff	[str]	Input the gff file
	
	-AAchange	[str] Input the SNP file that change the amino acid
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time
	
	2022.07.21 16:41 0.0.1
	
=cut

use strict;
use warnings;
use Getopt::Long;

my $analysis_path = "";
my ($snp, $ins, $del, $AAchange, $gff, $help) = ("", "", "", "", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'snp=s' => \$snp,
	'ins=s' => \$ins,
	'del=s' => \$del,
	'AAchange:s' => \$AAchange,
	'gff:s' => \$gff,
	'h|help:s' => \$help,
);

die `pod2text $0` if((!$analysis_path) or (!$snp) or (!$ins) or (!$del) or ($help));
die `pod2text $0` if(((!$gff) and ($AAchange)) or ((!$AAchange) and ($gff)));

my %SNPinGene;
if($AAchange){
	open(FILE, "<", $analysis_path.$AAchange);
	while(<FILE>){
		chomp;
		my @line = split/\t/;
		if($line[0] eq "Site"){
			next;
		}
		else{
			$SNPinGene{$line[0]} = join("\t", ($line[0], "SNP", @line[1..$#line]))."\n";
		}
	}
	close FILE;
}

my @out;
if($AAchange){
	push @out, join("\t", ("Site", "Mutation type", "Mutation", "Mutation in the gene", "Amino acid change", "Synonymous or non-syn"))."\n";
}
else{
	push @out, join("\t", ("Site", "Mutation type", "Mutation"))."\n";
}

open(FILE, "<", $analysis_path.$snp);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		if($AAchange){
			if($SNPinGene{$line[0]}){
				push @out, $SNPinGene{$line[0]};
			}
			else{
				push @out, join("\t", ($line[0], "SNP", $line[1], "Non-conding region", "No amino acid change", "NA"))."\n";
			}
		}
		else{
			push @out, join("\t", ($line[0], "SNP", $line[1]))."\n";
		}
	}
}
close FILE;

my %gene_hash;
if($gff){
	open(FILE, "<", $analysis_path.$gff);
	while(<FILE>){
		chomp;
		if($_ =~ m/#/ or $_=~/^\s*$/){
			next;
		}
		else{
			my @line = split/\t/,$_;
			my @gene = split/\;/,$line[$#line];
			if($line[2] eq "gene"){
				my $name = "FALSE";
				for my $i (@gene){
					if($i =~ m/Name=/){
						$name = "TRUE";
						$i =~ s/Name=//;
						@{$gene_hash{$i}} = ($line[3], $line[4], $line[6]);
						last;
					}
				}
				if($name eq "FALSE"){
					for my $i (@gene){
						if($i !~ m/Name=/ and $i =~ m/ID=/){
							$i =~ s/ID=//;
							@{$gene_hash{$i}} = ($line[3], $line[4], $line[6]);
							last;
						}
					}
				}
			}
		}
	}
	close FILE;
}

open(FILE, "<", $analysis_path.$del);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		my @delsites = split/\~/,$line[0];
		my @delseq = split/\:/,$line[1];
		if($AAchange){
			my $res = "FALSE";
			for my $gene (keys %gene_hash){
				if($delsites[0] >= @{$gene_hash{$gene}}[0] and $delsites[0] <= @{$gene_hash{$gene}}[1]){
					$res = "TRUE";
					my $delete_record;
					my ($delsite_start, $delsite_end);
					if(@{$gene_hash{$gene}}[2] eq "-"){
						my $gene_len = @{$gene_hash{$gene}}[1] - @{$gene_hash{$gene}}[0] + 1;
						$delsite_start = $gene_len - ($delsites[1] - @{$gene_hash{$gene}}[0] + 1) + 1;
						$delsite_end = $gene_len - ($delsites[0] - @{$gene_hash{$gene}}[0] + 1) + 1;
						my $delrc = reverse($delseq[1]);
						$delrc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
						$delete_record = "del:".$delsites[0]."-".$delseq[1]."-".$delsites[1]."(".$delsite_start."-".$delrc."-".$delsite_end.")";
					}
					else{
						$delsite_start = $delsites[0] - @{$gene_hash{$gene}}[0] + 1;
						$delsite_end = $delsites[1] - @{$gene_hash{$gene}}[0] + 1;
						$delete_record = "del:".$delsite_start."-".$delseq[1]."-".$delsite_end;
					}
					push @out, join("\t", ($line[0], "deletion", $delete_record, $gene, "-", "non-synonymous"))."\n";
					last;
				}
			}
			if($res eq "FALSE"){
				push @out, join("\t", ($line[0], "deletion", "del:".$delsites[0]."-".$delseq[1]."-".$delsites[1], "Non-conding region", "No amino acid change", "NA"))."\n";
			}
		}
		else{
			push @out, join("\t", ($line[0], "deletion", $line[1]))."\n";
		}
	}
}
close FILE;

open(FILE, "<", $analysis_path.$ins);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		my @insseq = split/\:/,$line[1];
		if($AAchange){
			my $res = "FALSE";
			for my $gene (keys %gene_hash){
				if($line[0] >= @{$gene_hash{$gene}}[0] and $line[0] <= @{$gene_hash{$gene}}[1]){
					$res = "TRUE";
					my $ins_record;
					my $ins_site;
					if(@{$gene_hash{$gene}}[2] eq "-"){
						my $gene_len = @{$gene_hash{$gene}}[1] - @{$gene_hash{$gene}}[0] + 1;
						$ins_site = $gene_len - ($line[0] - @{$gene_hash{$gene}}[0] + 1) + 1;
						my $insseq_rc = reverse($insseq[1]);
						$insseq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
						$ins_record = "ins:".$line[0]."-".$insseq[1]."(".$ins_site."-".$insseq_rc.")";
					}
					else{
						$ins_site = $line[0] - @{$gene_hash{$gene}}[0] + 1;
						$ins_record = "ins:".$ins_site."-".$insseq[1];
					}
					push @out, join("\t", ($line[0], "insertion", $ins_record, $gene, "-", "non-synonymous"))."\n";
					last;
				}
			}
			if($res eq "FALSE"){
				push @out, join("\t", ($line[0], "insertion", "ins:".$line[0]."-".$insseq[1], "Non-conding region", "No amino acid change", "NA"))."\n";
			}
		}
		else{
			push @out, join("\t", ($line[0], "insertion", $line[1]))."\n";
		}
	}
}
close FILE;
open(OUT, ">", $analysis_path."comparative-results.tsv");
print OUT @out;
close OUT;
