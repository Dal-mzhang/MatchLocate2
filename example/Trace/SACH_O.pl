#!/usr/bin/env perl
use warnings;
#Here we define 00h00m00.00s of one day as our reference time and mark the origin time as ZERO.
@ARGV == 1 || die"perl $0 dir\n";
$dir = $ARGV[0];
chomp($dir);

$hour = "00";
$min = "00";
$sec0 = "00";
$msec = "00";

@sac = glob "$dir/*";
open(SAC,"|sac> jk");
for($i=0;$i<@sac;$i++){
	chomp($sac[$i]);
	($jk,$year,$jday) = split(' ',`saclst nzyear nzjday f $sac[$i]`);chomp($jday);
	print SAC "rh $sac[$i]\n";
	print SAC "ch o gmt $year $jday $hour $min $sec0 $msec\n";
	print SAC "evaluate to tmp 0 - &1,o\n";
	print SAC "ch allt %tmp% iztype io\n";
	print SAC "ch o 0\n";
	print SAC "wh\n";
}
print SAC "q\n";
close(SAC);
unlink "jk";
