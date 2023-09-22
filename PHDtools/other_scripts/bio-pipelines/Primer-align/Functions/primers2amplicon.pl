#!/usr/bin/perl -w

=head1 Description
	Find the target sequence for a certain primers/probe set
Usage
	perl primers2amplicon.pl -WorkPath <path> -genomes <genome of target pathogen fasta> -primers <primers-probe fasta> -primerm6 <m6 file>
Parameters
	-WorkPath	[str]	Input the analysis path
	-genomes	[str]	Input the isolate genomes fasta
	-primers	[str]	Input the primers-probe fasta
	-primerm6	[str]	Input the m6 file of rapid primer alignment
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.18 21:21 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $analysis_path = "";
my ($isolates, $PCR_primers, $rapidm6, $help) = ("", "", "", "");

GetOptions(
	'WorkPath=s' => \$analysis_path,
	'genomes=s' => \$isolates,
	'primers=s' => \$PCR_primers,
	'primerm6=s' => \$rapidm6,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$isolates) or (!$PCR_primers) or (!$rapidm6) or ($help));

my (%direct, %bitscore, %info);
open(FILE, "<", $analysis_path.$rapidm6);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[8] < $line[9]){
		$direct{$line[0]}{"F"} += 1;
	}
	else{
		$direct{$line[0]}{"R"} += 1;
	}
	$bitscore{$line[1]} += $line[$#line];
	if(uc($line[0]) eq uc("forward")){
		@{$info{$line[1]}{"forward"}} = @line[8..9];
	}
	elsif(uc($line[0]) eq uc("reverse")){
		@{$info{$line[1]}{"reverse"}} = @line[8..9];
	}
}
close FILE;
my $primers = Bio::SeqIO -> new(-file => $analysis_path.$PCR_primers, -format => 'fasta');
my @pmlist;
while(my $primer_seq = $primers -> next_seq){
	my ($id, $seq_seq) = ($primer_seq -> id, $primer_seq -> seq);
	if(!$direct{$id}{"R"}){
		push @pmlist, join("\t", ($id, $seq_seq, "F"))."\n";
	}
	elsif(!$direct{$id}{"F"}){
		push @pmlist, join("\t", ($id, $seq_seq, "R"))."\n";
	}
	else{
		if($direct{$id}{"F"} > $direct{$id}{"R"}){
			push @pmlist, join("\t", ($id, $seq_seq, "F"))."\n";
		}
		else{
			push @pmlist, join("\t", ($id, $seq_seq, "R"))."\n";
		}
	}
}
open(OUT, ">", $analysis_path."pmlist.tsv");
print OUT @pmlist;
close OUT;
my $max_score = 0;
for my $id (sort{$bitscore{$b} <=> $bitscore{$a}} keys %bitscore){
	if($max_score > $bitscore{$id}){
		last;
	}
	else{
		$max_score = $bitscore{$id};
	}
}
my $amplicon;
for my $id (sort{$a cmp $b} keys %info){
	if($max_score == $bitscore{$id}){
		if(@{$info{$id}{"forward"}} and @{$info{$id}{"reverse"}}){
			if(@{$info{$id}{"forward"}}[0] < @{$info{$id}{"forward"}}[1] and @{$info{$id}{"reverse"}}[0] > @{$info{$id}{"reverse"}}[1]){
				my $iso_genome = Bio::SeqIO -> new(-file => $analysis_path.$isolates, -format => 'fasta');
				while(my $seq = $iso_genome -> next_seq){
					my ($seq_id, $seq_seq) = ($seq -> id, $seq -> seq);
					if($id eq $seq_id){
						$amplicon = substr($seq_seq, (@{$info{$id}{"forward"}}[0] -1 - 20), (@{$info{$id}{"reverse"}}[0] - @{$info{$id}{"forward"}}[0] + 1 + 40));
						last;
					}
				}
				last;
			}
		}
	}
}
open(OUT, ">", $analysis_path."AmpliconSeq.fasta");
print OUT ">target_sequence\n".$amplicon."\n";
close OUT;
