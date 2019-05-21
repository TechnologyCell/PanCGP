#!/usr/bin/perl
use strict;
#use Digest::MD5;
use Getopt::Long;

my ($dirpath);
&GetOptions (
# main options
 "dirpath:s"   => \$dirpath,       # keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )
);

newgenes();

sub newgenes {

    my @checksumlist;

    warn "#\nGROUPING CALLED\n";
    while ( <> ) {
        next if /^#/;
        chomp;
        my ($checksum,$totalgenes) = split /\./;
        print "description=$checksum, source=$totalgenes\n";
        my $id = $#checksumlist + 1;
        print (( $checksumlist[$id]->{checksum} , $checksumlist[$id]->{totalgenes} ) = ( $checksum , $totalgenes));
    }
    print "\nparse checksumlist end \n";

    open(my $FileWriter, ">", "$dirpath/newgenes.txt") or die "Could not open file $!";

    warn "#\nCALCULATING NEW GENES\n";
    foreach my $a ( 0 .. $#checksumlist ) {
        my %HITS;
        $checksumlist[$a]->{'newgenes'} = 0;
        print STDERR "# parsing blast reports (id $a)\n";
        foreach my $b ( 0 .. $a ) {
            open GUNZIP , "gunzip -c $dirpath/$a-$b.blastout.gz $dirpath/$b-$a.blastout.gz |";
            while (<GUNZIP>) {
                chomp;
                my ($q,$s) = split /\t/;
                $HITS{$q}{$s} = 1;
            }
            close GUNZIP;
        }
        open FSA , "$dirpath/$a.fsa" or die $!;
        while (<FSA>) {
            chomp;
            next unless /^>(.*)/;
            my $q = $1;
            my $hits_in_other_samples;
            foreach (keys %{$HITS{$q}}) {
                next unless /^(.*)\./; 
                if ( $checksumlist[$a]->{checksum} ne $1) {
                    $hits_in_other_samples = 1;
                    last;
                }
            }
            $checksumlist[$a]->{'newgenes'}++ if ! defined $hits_in_other_samples;
        }
        print $FileWriter $checksumlist[$a]->{'newgenes'}."\n";
    }
    close ($FileWriter);
}