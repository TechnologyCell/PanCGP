#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($dirpath);
&GetOptions (
 "dirpath:s"   => \$dirpath,       # keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )
);

coregenes();

sub coregenes {

    my @checksumlist;
    my $totallines = 0;

    while ( <> ) {
        $totallines++;
    }
    print "\nparse checksumlist end == $totallines \n";

    open (my $FileWriter, ">", "$dirpath/coregenes.txt") or die "Could not open file $!";

    foreach my $id (0 .. $totallines-1) {
        $checksumlist[$id]->{'coregenes'} = 0;

        open FileReader , "$dirpath/group_$id.dat";

        while (<FileReader>) {
            my ($group_id,$count,@members) = split /\t/;
            my %repr;
            foreach (@members) {
                    next unless /^(.*)\.(.*)/;
                    $repr{$1} = 1;
            }
            $checksumlist[$id]->{'coregenes'}++ if scalar ( keys %repr ) == $id + 1;
        }
        print "\n $checksumlist[$id]->{'coregenes'} \n";
        print $FileWriter $checksumlist[$id]->{'coregenes'}."\n";
        close FileReader;
    }
    close ($FileWriter);
}