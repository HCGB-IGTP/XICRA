#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0];
my $blocks = $ARGV[1];
my $ext = $ARGV[2];
my $folder = $ARGV[3];

if (scalar @ARGV == 0) {
        print "Usage:\nperl $0 file parts extension path\n\n";
        print "This script splits file in as many parts as stated...\n";
        exit();
}

my $size = -s $file; #To get only size
#print "+ Stats for file: $file\n";
#print "+ Chars: $size\n";
my $block = int($size/$blocks);

# Splits a file such a sam or whatever file that could be read for each line
open (FH, "<$file") or die "Could not open file $file";
#print "+ Splitting file $file into blocks of $block characters...\n";
my @files; my $j = 0;
while (1) {
   	my $chunk;
   	my @file_path = split("/", $file);
   	my @tmp = split (".".$ext, $file_path[-1]);
	my $file_name = $tmp[0];
   	my $block_file = $folder."/".$file_name."_part-".$j."_tmp.".$ext;
   	
	print $block_file."\n";
   	push (@files, $block_file);
   	open(OUT, ">$block_file") or die "Could not open destination file [DOMINO.pm:file_splitter]";
   	$j++;
   	if (!eof(FH)) { read(FH, $chunk,$block);  print OUT $chunk; } ## Print the amount of chars
   	if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken
   	close(OUT); last if eof(FH);
}
close(FH);
