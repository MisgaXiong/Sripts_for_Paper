#!/usr/bin/perl -w

=head1 Description
	Report the WGS quality
Usage
	perl covANDdep.pl -genome <genome of target pathogen> -DEPfile <depth file>
Parameters
	-genome	[str]	Input the name of forward fastq file
	-DEPfile	[str]	Input the name of reverse fastq file
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.13 16:25 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use POSIX;
use Getopt::Long;

my ($input_seq, $input_dep, $help) = ("", "", "");

GetOptions(
	'genome=s' => \$input_seq,
	'DEPfile=s' => \$input_dep,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$input_seq) or (!$input_dep) or ($help));

my $fasta = Bio::SeqIO -> new(
	-file => $input_seq,
	-format => 'fasta',
);

my (%seqs, %depths, %avgs);

while(my $seq = $fasta -> next_seq){
	my ($id, $seq_len) = ($seq -> id, $seq -> length);
	$seqs{$id} = $seq_len;
}

open(FILE, "<", $input_dep);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	$depths{$line[0]}{$line[1]} = $line[2];
	$avgs{"dep"}{$line[0]} += $line[2];
	if($line[2] >= 1){
		$avgs{"cov"}{$line[0]} += 1;
	}
}
close FILE;

my @out;
for my $id (sort{$a cmp $b} keys %seqs){
	my $avg_dep = $avgs{"dep"}{$id} / $seqs{$id};
	my $avg_cov = $avgs{"cov"}{$id} / $seqs{$id};
	my $window = 100;
	push @out, join("\t", ($id, $seqs{$id}."bp"));
	push @out, join("\t", ("Average depth", $avg_dep));
	push @out, join("\t", ("Coverage", $avg_cov));
	push @out, join("\t", ("Region", "Average depth"));
	if($seqs{$id} >= 200){
		my @tbl;
		my $len = ceil($seqs{$id}/$window);
		my $j = 1;
		my $win_dep = 0;
		for my $i (1..$seqs{$id}){
			if($j == $len and $i < $seqs{$id}){
				$win_dep += $depths{$id}{$i};
				$win_dep /= $len;
				push @tbl, join("\t", (($i - $len + 1), ($i), $win_dep));
				$win_dep = 0;
				$j = 1;
			}
			elsif($i == $seqs{$id}){
				$win_dep += $depths{$id}{$i};
				$win_dep /= $j;
				push @tbl, join("\t", (($i - $j + 1), ($i), $win_dep));
			}
			else{
				$win_dep += $depths{$id}{$i};
				$j += 1;
			}
		}
		push @out, join("\n", @tbl);
	}
	else{
		push @out, join("\t", (1, $seqs{$id}, $avg_dep));
	}
	push @out, "~" x 60;
}

print join("\n", @out)."\n"
