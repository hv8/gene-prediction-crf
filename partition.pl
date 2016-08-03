#!/usr/bin/env perl

use strict;

open FILELIST, "<random_filelist.txt";

print "Partitioning... ";
my @list_of_files; 
chomp(@list_of_files = <FILELIST>);

close FILELIST;

my $total_files = $#list_of_files + 1;
my $train_total = int($total_files * 0.7);

open OUTPUT, ">train.txt";
my $count = 0;
while ($count < $train_total)
{
	open FILE, "<$list_of_files[$count]";
	chomp(my @files = <FILE>);
	close FILE;
		
	foreach my $line (@files)
	{
		print OUTPUT "$line \n";
	}
	print OUTPUT "\n";
	
	$count++;	
}

close OUTPUT;

open OUTPUT, ">test.txt";
while ($count < $total_files)
{
	open FILE, "<$list_of_files[$count]";
	chomp(my @files = <FILE>);
	close FILE;
		
	foreach my $line (@files)
	{
		print OUTPUT "$line \n";
	}
	print OUTPUT "\n";
	
	$count++;	
}

close OUTPUT;