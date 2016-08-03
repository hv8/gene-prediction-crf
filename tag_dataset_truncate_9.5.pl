#!/usr/bin/perl

########################
##  Global Variables  ##
########################

$base_index;
@bases;
$base_pairs;
$C_count;
@fsm_output;

@histo_output;
@histo_output_a;
@histo_output_t;
@histo_output_c;
@histo_output_g;
@histo_output_gt;
@histo_output_ag;
@histo_output_atg;
@histo_output_taa;
@histo_output_tga;
@histo_output_tag;
@histo_output_a_f;
@histo_output_t_f;
@histo_output_c_f;
@histo_output_g_f;
@histo_output_gt_f;
@histo_output_ag_f;
@histo_output_atg_f;
@histo_output_taa_f;
@histo_output_tga_f;
@histo_output_tag_f;

@tags_output;
@join_data;

open FILELIST, "<filelist.txt";

my @list_of_files; 
chomp(@list_of_files = <FILELIST>);
close FILELIST;

open STATFILE, ">>counts2.txt";
foreach $filename (@list_of_files)
{
    $outfilename = "dataset_".$filename;

    open(DATASET, "+<", "$filename");
    open(OUTFILE, ">$outfilename");
    &gene_concat;
    &format_gene_sequence;
    &data_label;
    
    close(OUTFILE);
    close(DATASET);
    
    print STATFILE "$filename \t $base_pairs \t $C_count\n";
    print "$filename $base_pairs $C_count\n";
        
    $base_index = 0;
    @bases = ();
    $base_pairs = 0;
    @fsm_output = ();
    $C_count = 0;

    @histo_output = ();
    @histo_output_a = ();
    @histo_output_t = ();
    @histo_output_c = ();
    @histo_output_g = ();
    @histo_output_gt = ();
    @histo_output_ag = ();
    @histo_output_atg = ();
    @histo_output_taa = ();
    @histo_output_tga = ();
    @histo_output_tag = ();
    @histo_output_a_f = ();
    @histo_output_t_f = ();
    @histo_output_c_f = ();
    @histo_output_g_f = ();
    @histo_output_gt_f = ();
    @histo_output_ag_f = ();
    @histo_output_atg_f = ();
    @histo_output_taa_f = ();
    @histo_output_tga_f = ();
    @histo_output_tag_f = ();

    @tags_output = ();
    @join_data = ();
}
close STATFILE;

##########################
##  Gene Concatenation  ##
##########################
sub gene_concat
{
    $gene_sequence="";
    $line;
    while($line = <DATASET>)
    {	
            chomp($line);
            @base_array = split("", $line);
            if(@base_array[0] =~ /[GTCA]/)
            {
                    $gene_sequence .= $line
            }
            elsif(@base_array[0] =~ /[j]/)
            {
                    # parse join data
            &parse_join_data;
            }
    }
}


##############################
## Formatting Gene Sequence ##
##############################
sub format_gene_sequence
{
    @bases=split("",$gene_sequence);
    $base_pairs=@bases;

    for($i=0; $i<$base_pairs; $i++)
    {
            if(@bases[$i] =~ /[^GTCA]/)
            {
                    splice @bases, $i, 1;
                    $base_pairs--;
                    $i--;
            }

    }

    $base_pairs = @bases;
    printf "\nGene has a total of %d bases.\n\n", $base_pairs;
}
#########################
##  Dataset Labeling   ##
######################### 
sub data_label
{

    ## truncate bases ##
    &truncate_introns;

    $base_pairs = @bases;
    printf "\nAfter truncation: Gene has a total of %d bases.\n\n", $base_pairs;

    ## reset ##
    @tags_output=();
    $base_pairs = @bases;

    ## tag each base ##
    &gene_tagging;
        
    #&fsm_feat;
    &histogram_feat;

    $tags_len=@tags_output;
    printf "\nTotal Tags length = %d\n", $tags_len;

    $histo_g_len = @histo_output_g;
    $histo_gt_len = @histo_output_gt;
    $histo_tag_len = @histo_output_tag;
    printf "Histogram G Counts Length = %d\n", $histo_g_len;
    printf "Histogram GT Counts Length = %d\n", $histo_gt_len;
    printf "Histogram TAG Counts Length = %d\n", $histo_tag_len;

    #writing output arrays to file
    for($i=0; $i < $base_pairs; $i++)
    {
        printf OUTFILE "@bases[$i] %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d @tags_output[$i]\n", @histo_output_a[$i], @histo_output_t[$i], @histo_output_c[$i], @histo_output_g[$i], @histo_output_gt[$i], @histo_output_ag[$i], @histo_output_atg[$i], @histo_output_taa[$i], @histo_output_tga[$i], @histo_output_tag[$i], @histo_output_a_f[$i], @histo_output_t_f[$i], @histo_output_c_f[$i], @histo_output_g_f[$i], @histo_output_gt_f[$i], @histo_output_ag_f[$i], @histo_output_atg_f[$i], @histo_output_taa_f[$i], @histo_output_tga_f[$i], @histo_output_tag_f[$i];
    }
    printf OUTFILE "\n";

    print "\nDone labelling the bases of this gene.\n\n";
}

####################################
##   Subroutine for Gene Tagging  ##
####################################

sub gene_tagging
{
    $num_join_pairs = @join_data;
    $C_count = 0;
    printf "Join Indexes:\n";
    for($i=0; $i < $num_join_pairs; $i++)
    {
        printf "%d - %d\n", $i, @join_data[$i];
    }
    
	for($master_idx=0, $join_data_index=0; $master_idx < $base_pairs && $join_data_index <= $num_join_pairs; )
	{		
		#beginning of the sequence, label bases as N
		if($join_data_index==0)
                {
			for($j=$master_idx; $j < @join_data[$join_data_index]; $j++)
			{
				#print OUTFILE2 "@bases[$j] N\n";
                                push(@tags_output, "N");
                        }
			$master_idx = @join_data[$join_data_index];
			$join_data_index++;
		}
        
                #end of the sequence, label bases as N
                elsif($join_data_index==$num_join_pairs)
                {
                    for($j=$master_idx; $j < $base_pairs; $j++)
                    {
                        #print OUTFILE2 "@bases[$j] N\n";
                        push(@tags_output, "N");
                    }
                    $master_idx = $base_pairs;
                }
                
                #middle of the sequence, label bases as C or N
                else
                {
                    if($join_data_index%2 != 0) #label bases as C
                    {
                        printf "Start of coding: %d\n", $master_idx;
                        $tag_data_len = @tags_output;
                        printf "Length of Tag data so far: %d\n", $tag_data_len;
                        for($j=$master_idx; $j<=@join_data[$join_data_index]; $j++)
                        {
                            #print OUTFILE2 "@bases[$j] C\n";
                            push(@tags_output, "C");
                            $C_count++;
                        }
                        $master_idx=@join_data[$join_data_index]+1;
                        $join_data_index++;
                     }
                     
                     elsif($join_data_index%2 ==0) #label bases as N
                     {
                        for($j=$master_idx; $j<@join_data[$join_data_index]; $j++)
                        {
                            #print OUTFILE2 "@bases[$j] N\n";
                            push(@tags_output, "N");
                        }
                        $master_idx=$master_idx=@join_data[$join_data_index];
                        $join_data_index++;
                     }
                }
	}

}

########################################
##  Subroutine for Histogram Feature  ##
########################################

sub histogram_feat
{
    $gt_count=0;
    $ag_count=0;
    $atg_count=0;
    $taa_count=0;
    $tga_count=0;
    $tag_count=0;
    
    $g_count=0;
    $c_count=0;
    $t_count=0;
    $a_count=0;    
    
    $gt_f_count=0;
    $ag_f_count=0;
    $atg_f_count=0;
    $taa_f_count=0;
    $tga_f_count=0;
    $tag_f_count=0;
    
    $g_f_count=0;
    $c_f_count=0;
    $t_f_count=0;
    $a_f_count=0;
    
	$history_size=40;
	$future_size=20;
	
	$sequence_length = $#bases + 1;
    
    ## singleton and conjunction counts for previous 40 bases ##
    $index=0;
    #$history_window_count=int($base_pairs/$history_size);
    #printf "No. of blocks = %d\n", $history_window_count;
    for($i=0; $i<$base_pairs; $i++)
    {
        my $j = 1;
        while (($j <= $history_size) && ($j <= $i))
        {
            $index = $i - $j;
            
            #check for singletons
            if(@bases[$index] eq 'G')
            {
                $g_count++;
            }
            elsif(@bases[$index] eq 'C'){
                $c_count++;
            }
            elsif(@bases[$index] eq 'T'){
                $t_count++;
            }
            elsif(@bases[$index] eq 'A'){
                $a_count++;
            }
            
            
            #check for conjunctions of size 2
            if(@bases[$index-1] eq 'G' && @bases[$index] eq 'T')
            {
                $gt_count++;
            }
            elsif(@bases[$index-1] eq 'A' && @bases[$index] eq 'G')
            {
                $ag_count++;
            }
            
            
            #check for conjunctions of size 3
            if(@bases[$index-2] eq 'A' && @bases[$index-1] eq 'T' && @bases[$index] eq 'G')
            {
                $atg_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'A' && @bases[$index] eq 'A')
            {
                $taa_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'G' && @bases[$index] eq 'A')
            {
                $tga_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'A' && @bases[$index] eq 'G')
            {
                $tag_count++;
            }
            
            $j++;
        }
     
        my $j = 1;
        while (($j <= $future_size) && ($i + $j < $base_pairs))
        {
            $index = $i + $j;
            
            #check for singletons
            if(@bases[$index] eq 'G')
            {
                $g_f_count++;
            }
            elsif(@bases[$index] eq 'C'){
                $c_f_count++;
            }
            elsif(@bases[$index] eq 'T'){
                $t_f_count++;
            }
            elsif(@bases[$index] eq 'A'){
                $a_f_count++;
            }
            
            
            #check for conjunctions of size 2
            if(@bases[$index-1] eq 'G' && @bases[$index] eq 'T')
            {
                $gt_f_count++;
            }
            elsif(@bases[$index-1] eq 'A' && @bases[$index] eq 'G')
            {
                $ag_f_count++;
            }
            
            
            #check for conjunctions of size 3
            if(@bases[$index-2] eq 'A' && @bases[$index-1] eq 'T' && @bases[$index] eq 'G')
            {
                $atg_f_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'A' && @bases[$index] eq 'A')
            {
                $taa_f_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'G' && @bases[$index] eq 'A')
            {
                $tga_f_count++;
            }
            elsif(@bases[$index-2] eq 'T' && @bases[$index-1] eq 'A' && @bases[$index] eq 'G')
            {
                $tag_f_count++;
            }
            
            $j++;
        }
      
       
            push(@histo_output_a, $a_count);
            push(@histo_output_t, $t_count);
            push(@histo_output_c, $c_count);
            push(@histo_output_g, $g_count);
            push(@histo_output_gt, $gt_count);
            push(@histo_output_atg, $atg_count);
            push(@histo_output_taa, $taa_count);
            push(@histo_output_tga, $tga_count);
            push(@histo_output_tag, $tag_count);
            
            push(@histo_output_a_f, $a_f_count);
            push(@histo_output_t_f, $t_f_count);
            push(@histo_output_c_f, $c_f_count);
            push(@histo_output_g_f, $g_f_count);
            push(@histo_output_gt_f, $gt_f_count);
            push(@histo_output_atg_f, $atg_f_count);
            push(@histo_output_taa_f, $taa_f_count);
            push(@histo_output_tga_f, $tga_f_count);
            push(@histo_output_tag_f, $tag_f_count);     
        
        #reset the counts for the next block
        $gt_count=0;
        $ag_count=0;
        $atg_count=0;
        $taa_count=0;
        $tga_count=0;
        $tag_count=0;
        
        $g_count=0;
        $c_count=0;
        $t_count=0;
        $a_count=0;
    
    $gt_f_count=0;
    $ag_f_count=0;
    $atg_f_count=0;
    $taa_f_count=0;
    $tga_f_count=0;
    $tag_f_count=0;
    
    $g_f_count=0;
    $c_f_count=0;
    $t_f_count=0;
    $a_f_count=0;
    }
}



###################################
## Subroutine to parse join data ##
###################################

sub parse_join_data()
{
    $join_data_str = $line;
    
    @join_array=split(",", $join_data_str);
    $join_array_len=@join_array;
    
    for($i=0; $i<$join_array_len; $i++)
    {        
        @join_elements=split(/\.\./, @join_array[$i]);
        
        $join_elem_len = @join_elements;
        push(@join_data, @join_elements[0]);
        push(@join_data, @join_elements[1]);
    }
    
    $extract_len = length(@join_data[0])-5;
    $extracted_substr = substr(@join_data[0], 5, $extract_len);
    @join_data[0] = $extracted_substr;
    
    $last_index = @join_data;
    $last_index-=1;
    $extract_len = length(@join_data[$last_index]);
    $extract_len-=1;
    $extracted_substr = substr(@join_data[$last_index], 0, $extract_len);
    @join_data[$last_index] = $extracted_substr;
    
    $join_data_length=@join_data;
    for($i=0; $i<$join_data_length; $i++)
    {
        @join_data[$i]=@join_data[$i]-1;
    }
}



#######################################
## Subroutine for Truncating Introns ##
#######################################

sub truncate_introns
{    
    $num_join_pairs = @join_data;
    
    printf "Before truncation: join indexes:\n";
    for($i=0; $i < $num_join_pairs; $i++)
    {
        printf "%d - %d\n", $i, @join_data[$i];
    }
    
    $max_n_length=200;
    
    my $n_length = $join_data[0];
    my $portion_to_cut; 
       
    if ($n_length > 200)
    {
        my $second_portion =  $join_data[0] - 100;
        $portion_to_cut = $second_portion - 100;
        splice @bases, 100, $portion_to_cut; 
        
    
        for ($i=0; $i < $num_join_pairs; $i++)
        {
            $join_data[$i] = $join_data[$i] - $portion_to_cut;
        }
        $base_pairs = $base_pairs - $portion_to_cut;    
    }


    my $i = 1;
    while (($n_length < 200) && ($i < $num_join_pairs))
    {
        $n_length = $join_data[$i+1] - $join_data[$i] - 1;
        if ($n_length < 200)
        {
            $i+=2;
        }
    }
    for ( ; $i < $num_join_pairs; $i+=2)
    {
        $n_length = $join_data[$i+1] - $join_data[$i] - 1;
        if ($n_length > 200)
        {
            $second_portion =  $join_data[$i+1] - 100;
            $portion_to_cut = $second_portion - 100 - $join_data[$i] - 1;
            splice @bases, $join_data[$i] + 100, $portion_to_cut;
                    
            my $j = $i+1;
            for ( ; $j < $num_join_pairs; $j++)
            {
                #print "$j \n";
                #print "$join_data[$j] \n";
                $join_data[$j] = $join_data[$j] - $portion_to_cut;
                #print "$join_data[$j] \n";
            }
            $base_pairs = $base_pairs - $portion_to_cut;   
        } 
    }        

    $i = $num_join_pairs - 1;
    print "i = $i \n";
    $n_length = $base_pairs - $join_data[$i] - 1;
    if ($n_length > 200)
    {
            $second_portion =  $base_pairs - 100;
            $portion_to_cut = $second_portion - 100 - $join_data[$i] - 1;
            splice @bases, $join_data[$i] + 100, $portion_to_cut;
            $base_pairs = $base_pairs - $portion_to_cut;        
    }
}



