#!/usr/bin/perl -w

=head1 Description
	Delete the PHDtools results after a period of time
Usage
	perl deletdata.pl -nowdate <current data year-month-day-hour-min-second>
Parameters
	-nowdate	[str]	Input the current date
	-h/-help	[str]	print help
Auther
	Dongyan Xiong
Edit Time
	2022.08.30 22:42 0.0.1
=cut

use strict;
use warnings;
use Getopt::Long;

my ($current_date, $help) = ("", "");
GetOptions(
	'nowdate=s' => \$current_date,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$current_date) or ($help));

my $time = $current_date;
my @times = split/\-/,$time;
my $storage_path = "/home/dell/data/Web-Server/Result-store/";
my @folders = `ls -R $storage_path | grep \.\/`;
for my $folder (@folders){
	$folder =~ s/\n//;
	$folder =~ s/\://;
	my @parent = split/\//,$storage_path;
	my @son = split/\//, $folder;
	if($#son - $#parent == 2){
		my $target_folder = $son[$#son];
		my @comp = split/\-/,$target_folder;
		my $res = "TRUE";
		for my $var (@comp){
			if($var =~ /[^0-9]+/){
				$res = "FALSE";
			}
		}
		if($res eq "TRUE"){
			my ($year, $month, $date) = (($times[0] - $comp[0]), ($times[1] - $comp[1]), ($times[2] - $comp[2]));
			if($year >= 1 or $month >= 1 or $date >= 5){
				my $command = "rm -rf ".$folder;
				system($command);
			}
		}
	}
}
