#!/usr/local/bin/perl

$filename = @ARGV[0];
$first    = @ARGV[1];
$last     = @ARGV[2];

open ( OUT, ">$first-$last$filename" ) || die "Errore su $filename";

$i = $first; 
for($i; $i<$last; $i++){
    print OUT ("/nfs/farm/babar/AWG10/LeptBc/BtauNuSemiExcl/production/tautau/tautau-$i.root\n");
}
print OUT ("/nfs/farm/babar/AWG10/LeptBc/BtauNuSemiExcl/production/tautau/tautau-$last.root");
close(OUT);

