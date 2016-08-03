#!/usr/bin/env perl

use strict;

open RESULT, "<results.txt";
my @results;
my $total = 0;
while (my $line = <RESULT>)
{
	chomp $line;
	$line =~ s/\0//ig;

	unless ($line =~ /^\s*$/) { 
		push(@results, $line);
		$total++;
	}
	else
	{
		push(@results, " ");
	}

}
close RESULT;

my $TP = 0;
my $FN = 0;
my $TN = 0;
my $FP = 0;
my $total_C = 0;
my $total_N = 0;
my $totalgenes = 0;
		
foreach my $line (@results)
{
	my $line_length = length($line);
	if ($line ne " ")
	{
		my @result_line = split/\s+/, $line;
		my ($true_value, $predicted) = ($result_line[21], $result_line[22]);
			
		if ($true_value eq 'C')
		{
			if ($predicted eq 'C')
			{
				$TP++;
			}
			else
			{
				$FN++;
			}
			$total_C++;
		}
		if ($true_value eq 'N')
		{
			if ($predicted eq 'N')
			{
				$TN++;
			}
			else
			{
				$FP++;
			}
			$total_N++;			
		}		
	}
}

print "TP = $TP \n";
print "FN = $FN \n";
print "TN = $TN \n";
print "FP = $FP \n";
print "Total C = $total_C \n";
print "Total N = $total_N \n";

my $sensitivity = $TP / ($TP + $FN);
my $specificity = $TN / ($TN + $FP);
my $F1 = (2*$sensitivity* $specificity) / ($sensitivity + $specificity);

#print "Total Words = $totalgenes \n";
print "Total = $total \n";
print "Sensitivity = $sensitivity \n";
print "Specificity = $specificity \n";
print "F1 = $F1 \n";

open OUTPUT, ">>evaluate.txt";
print OUTPUT "TP = $TP \n";
print OUTPUT "FN = $FN \n";
print OUTPUT "TN = $TN \n";
print OUTPUT "FP = $FP \n";
print OUTPUT "Total C = $total_C \n";
print OUTPUT "Total N = $total_N \n";

#print "Total Words = $totalgenes \n";
print OUTPUT "Total = $total \n";
print OUTPUT "Sensitivity = $sensitivity \n";
print OUTPUT "Specificity = $specificity \n";
print OUTPUT "F1 = $F1 \n";
print OUTPUT "\n";
close OUTPUT;
