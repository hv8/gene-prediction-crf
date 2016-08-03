#!/usr/bin/env perl

use strict;

my @list_of_files;

my $count = 1;

while ($count < 10)
{
	push(@list_of_files, "dataset_00$count.txt");
	$count++;
}
while ($count < 50 + 1)
{
	push(@list_of_files, "dataset_0$count.txt");
	$count++;
}

my $i = $#list_of_files;
while (--$i)
{
	my $j = int rand ($i + 1);
	@list_of_files[$i, $j] = @list_of_files[$j, $i];
}

open OUTPUT, ">random_filelist.txt";
foreach my $file(@list_of_files)
{
	print "$file \n";
	print OUTPUT "$file \n";
}
close OUTPUT;