#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($dirpath);
&GetOptions (
# main options
 "dirpath:s"   => \$dirpath,       # keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )
);

newgenefamilies();

sub newgenefamilies {

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

    open(my $FileWriter, ">", "$dirpath/newgenefamilies.txt") or die "Could not open file $!";

    warn "#\nCALCULATING NEW GENES\n";
    foreach my $id (0 .. $#checksumlist) {
        $checksumlist[$id]->{'newgenefamilies'} = 0;
        open GRP , "$dirpath/group_$id.dat";
        while (<GRP>) {
            my ($group_id,$count,@members) = split /\t/;
            my $is_new_family = 1;
            foreach (@members) {
                next unless /^(.*)\.(.*)/;
                $is_new_family = 0 if $1 ne $checksumlist[$id]->{checksum};
            }
            $checksumlist[$id]->{'newgenefamilies'} += $is_new_family;
        }
        print $FileWriter $checksumlist[$id]->{'newgenefamilies'}."\n";
        close GRP;
    }
    close ($FileWriter);
}