#!/usr/bin/env perl

use strict;

open OUTPUT, ">datafilelist.txt";

my $count = 1;

while ($count < 10)
{
	print OUTPUT "dataset_00$count.txt\n";
	$count++;
}
while ($count < 50 + 1)
{
	print OUTPUT "dataset_0$count.txt\n";
	$count++;
}

close OUTPUT;