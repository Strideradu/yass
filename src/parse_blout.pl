#!/usr/bin/perl
use strict;
use Bio::SearchIO;

($#ARGV == 1) or die "Usage: $ARGV[0] blastoutputfile > inputfile\n";

my $searchio = new Bio::SearchIO(-format => 'blast',
				 -file   => $ARGV[0]);
while( my $result = $searchio->next_result ) {
    while( my $hit = $result->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
	    print "*(".$hsp->start('query')."-".$hsp->end('query').")-(".$hsp->start('subject')."-".$hsp->end('subject').") Ev: ".$hsp->evalue." s: ".$hsp->length('query')."/".$hsp->length('subject')."\n";
	}
    }
}
