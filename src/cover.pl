#!/usr/bin/perl -w 
use strict;
$| = 1;

# max two args
sub max2($$) {
    my ($x, $y) = @_;
    return (($x > $y) ? $x : $y);
}

# min two args
sub min2($$) {
    my ($x, $y) = @_;
    return (($x < $y) ? $x : $y);
}

# max several args
sub max($) {
    my $m = -0x80000000;
    while(scalar(@_)){
	$m = &max2($m,shift @_);
    };
    return $m;
}

# min several args
sub min($) {
    my $m = +0x7fffffff;
    while(scalar(@_)){
	$m = &min2($m,shift @_);
    };
    return $m;
}

# line to parse, output table to store
sub parseline($$){
    my ($line,$ptarray) = @_;
    m/.*\t.*\t.*\t.*\t.*\t.*\t([0-9]+)\t([0-9]+)\t([0-9]+)\t([0-9]+)\t.*\t.*/;
    @{$ptarray}[0] = &min($1,$2);
    @{$ptarray}[1] = &max($1,$2);
    @{$ptarray}[2] = &min($3,$4);
    @{$ptarray}[3] = &max($3,$4);
    #print " ".@{$ptarray}[0]."\t".@{$ptarray}[1]."\t".@{$ptarray}[2]."\t".@{$ptarray}[3]."\n"
}

# overlap between two 
sub overlap($$){
    my ($ptarray1,$ptarray2) = @_;    
    my  $res =  !( 
		(@{$ptarray1}[0] > @{$ptarray2}[1] || @{$ptarray2}[0] > @{$ptarray1}[1]) ||
		(@{$ptarray1}[2] > @{$ptarray2}[3] || @{$ptarray2}[2] > @{$ptarray1}[3]))
	;
    return $res;
}

# length 
sub len($){
    my ($line) = @_;
     m/.*\t.*\t.*\t([0-9]+)\t.*\t.*\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t.*\t.*$/;
    return $1;
}

# evalue
sub evalue($){
    my ($line) = @_;
    m/.*\t.*\t.*\t[0-9]+\t.*\t.*\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t(.*)\t.*$/;
    return $1;
}

# percent
sub percent($){
    my ($line) = @_;
    m/.*\t.*\t([0-9]*\.[0-9]*)\t.*\t.*\t.*\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t.*\t.*/;
    return $1;
}


# comparison procedure
sub difference($$$) {
    my ($evalue,$file1,$file2) = @_;
    open(IN1, $file1) or die("* cannot open \"$file1\"");
    
    my $r_total     = 0;
    my $r_uncovered = 0;
    my $r_len       = 0;
    my $r_total_len = 0;
# main loop 
    while(<IN1>){
	if (m/^\#.*$/ || m/^$/ ){
	    #first header line 
	}else{
	    #data line
	    if (&evalue($_) <= $evalue){
		my $len = &len($_);
		$r_total_len += $len;
		$r_total++;
		my @pos1 = (0,0,0,0);
		parseline($_,\@pos1);
		#print " ".$pos1[0]."\t".$pos1[1]."\t".$pos1[2]."\t".$pos1[3]."\n";
		open(IN2, $file2) or die("* cannot open \"$file2\"");	
		my $over = 0;
		while(<IN2>){	
		    if (m/^\#.*$/ || m/^$/ ){
			#first header line 
		    }else{
			#data line
			if (1){
			    my @pos2 = (0,0,0,0);
			    parseline($_,\@pos2);			
			    if (overlap(\@pos1, \@pos2 ) == 1) {
				$over = 1; 
				#print "#".$pos2[0]."\t".$pos2[1]."\t".$pos2[2]."\t".$pos2[3]."\n";
				goto e;
			    }
			}
		    }
		}
	      e:
		close(IN2);
		if ($over == 0){
		    $r_uncovered++;
		    $r_len += $len;
		}
	    }
	}
    }
    close(IN1);

# output number of repeats uncovered and number of repeats 
    print " missing repeats of \"".$file2."\" compared to reference set \"".$file1."\" with evalue < ".$evalue.")\n";
    print " number of uncovered : ".$r_uncovered."/".$r_total."\t (length : ".$r_len."/".$r_total_len.")\n";
}




#--------------------------------------------#
# MAIN

# parameters : 
# len1 percent1 len2 percent2 file1 file2
($#ARGV + 1 ==  2)  or die("* two blast output files are needed\n");

for (my $evalue  = 1e-9 ; $evalue <= 10.0 ; $evalue *= 10) {
    print("\n");
    &difference($evalue,$ARGV[0],$ARGV[1]);
    &difference($evalue,$ARGV[1],$ARGV[0]);
}



