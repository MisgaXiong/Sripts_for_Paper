#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $analysis_path = "";
my ($PCR_primers, $help) = ("", "");

GetOptions(
	'WorkPath=s' => \$analysis_path,
	'primers=s' => \$PCR_primers,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$PCR_primers) or ($help));

my $primers = Bio::SeqIO -> new(-file => $analysis_path.$PCR_primers, -format => 'fasta');
my $total = 0;
my @nPCR = ("forward", "reverse");
my @qPCR = ("forward", "reverse", "probe");
my %primerseqs;
while(my $seq = $primers -> next_seq){
	my ($id, $seq_seq) = ($seq -> id, $seq -> seq);
	$primerseqs{$id} = $seq_seq;
	$total += 1;
}
my $res = "TRUE";
if($total > 3 or $total < 2){
	$res = "FALSE";
}
elsif($total >= 2 and $total <= 3){
	if($total == 2){
		for my $id (@nPCR){
			if(!$primerseqs{$id}){
				$res = "FALSE";
				last;
			}
		}
	}
	elsif($total == 3){
		for my $id (@qPCR){
			if(!$primerseqs{$id}){
				$res = "FALSE";
				last;
			}
		}
	}
}

print $res;

