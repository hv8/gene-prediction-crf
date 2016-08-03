#!/usr/bin/env perl

use strict;

open OUTPUT, ">filelist.txt";

my $count = 1;

while ($count < 10)
{
	print OUTPUT "00$count.txt\n";
	$count++;
}

while ($count < 50 + 1)
{
	print OUTPUT "0$count.txt\n";
	$count++;
}

close OUTPUT;