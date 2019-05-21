#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($dirpath);
my ($rank);
&GetOptions (
# main options
 "dirpath:s"   => \$dirpath,       # keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )
 "rank:s"  => \$rank,
);

grouping();

sub grouping {
    my @checksumlist;
    #warn "#\nGROUPING CALLED\n";
    while ( <> ) {
        next if /^#/;
        chomp;
        my ($checksum,$totalgenes) = split /\./;
        print "description = $checksum, source = $totalgenes\n";
        my $id = $#checksumlist + 1;
        print (( $checksumlist[$id]->{checksum} , $checksumlist[$id]->{totalgenes} ) = ( $checksum , $totalgenes));
    }
    print "\nparse checksumlist end \n";

        open GRP , "| perl /home/ahsan/NetBeansProjects/mpjtest/src/mpjtest/group > $dirpath/group_$rank.dat";
        foreach my $x (0 .. $rank) {
            foreach my $y (0 .. $rank) {
                foreach my $z ($x , $y) {
                    foreach my $n ( 1 .. $checksumlist[$z]->{'totalgenes'} ) {
                        print GRP $checksumlist[$z]->{checksum}.".$n\n";
                    }
                }
                open BLASTOUT , "gunzip -c $dirpath/$x-$y.blastout.gz |";
                while (<BLASTOUT>) {
                    print GRP $_;
                }
                close BLASTOUT;
            }
        }
        close GRP;
}

