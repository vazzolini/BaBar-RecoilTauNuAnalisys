#!/usr/local/bin/perl -w


# argomenti : string ( tipo file root: generic cocktail etc... )
#           : maxi (max num di root della chain)
#           : filetype (nome chain di output)

# per esempio $string = "/nfs/farm/babar/AWG8/ISL/041302/root/genbch/d0-sp4run1-gen-bb-8000-";
# per esempio $filetype = "genbchchain_d0_run1_";

$filename=$ARGV[0];

$j=0;

open(FILIN,"$filename") || die "errore in $filename";
while (<FILIN>) {
( $burp, $rootto, $dummy )=split(/\: /,$_);
#print ("$burp\n");
print ("$rootto\n");
#print ("$dummy\n");
$j += $rootto;
}
close(FILIN);
print ("$j\n");










