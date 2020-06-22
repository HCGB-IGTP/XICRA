#!/usr/bin/perl
use strict;
use warnings;

##############################################################################
##	Jose Fco. Sanchez Herrero, 21/05/2019 jfbioinformatics@gmail.com		##
##############################################################################

my $ids = $ARGV[0];
my $file = $ARGV[1];
my $output = $ARGV[2];

if (!@ARGV) {
	print "############################\n";
	print "\tHELP Message:\n";
	print "############################\n";
	print "\nThis script retrieves the ids provided in a txt file from a fasta file (protein or nucleotide).\n";
	print "Usage:\nperl $0 ids.txt sequences.fasta output_file.fasta\n\n";
	exit();
}

## get ids
open (ID,"$ids");
my @IDS = (<ID>);
close (ID);
chomp @IDS;

#read proteins
my %hash;
open (OUT, ">$output");
open(FILE, $file) || die "Could not open the $file ...\n";
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {
	#next if /^#/ || /^\s*$/;
	chomp;
	my ($titleline, $sequence) = split(/\n/,$_,2);
	next unless ($sequence && $titleline);
	my @id = split(" ", $titleline); 
	chop $sequence;
	$sequence =~ s/\n//g;
	if (grep /$id[0]$/, @IDS) {
		print OUT ">".$titleline."\n".$sequence."\n";
	}
}
close(FILE); $/ = "\n";
close (OUT);

