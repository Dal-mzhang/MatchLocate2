#!/usr/bin/env perl
use strict;
use warnings;
#Usually, only direct waves (e.g., Pg, Sg) are used in a regional area.
my $dir = ".";
my @files = `ls -d $dir/2012*`;
#mutiple phase
my $mult=1;

#P wave
#if ts-tp is larger than $tdif, phase_p would be used.
#Otherwise, the largest phase S would be dominated in the P window.
#$tdif depends on your correlation window length.
my $phase_p = "P,p";
my $tdif = "4"; 

#S wave
my $phase_s = "S,s"; 
my $distmax = "200";

#my $model = "JMA"; If you have model JMA2001.
my $model = "prem";

my $INPUT = "INPUT";
if(-e $INPUT){}else{`mkdir $INPUT`;}

my $jk="jk";
open(SAC,"|sac >$jk");
for(my $i=0;$i<@files;$i++){
	chomp($files[$i]);
	my ($jk,$FILE) = split("\/",$files[$i]);
	my $NAME = "$INPUT/$FILE";
	my @sta = `ls $files[$i]/*`;
	open(FL,">$NAME");
foreach $_(@sta){
	chomp($_);
	my ($jk0,$gcarc,$dist,$kstnm,$kcmpnm,$evdp,$mag) = split(" ",`saclst gcarc dist kstnm kcmpnm evdp mag f $_`);chomp($mag);

	if($dist <=$distmax){
		print"$_\n";

		my ($tp,$dt_dgc_p,$dt_dh_p) = &taup($model,$evdp,$dist,$phase_p);
		my ($ts,$dt_dgc_s,$dt_dh_s) = &taup($model,$evdp,$dist,$phase_s);
		if($tp < 0 || $ts < 0){print STDERR "stop!! Wrong with the arrival time!\n"; exit;}
			print SAC "rh $_\n";
			print SAC "ch t1 -12345\nch t2 -12345\nch user0 -12345\n";
			print SAC "wh\n";
		#use multiple phases to constrain the location (e.g., P and S, P and Pp, at al.).
		if((($ts - $tp) > $tdif) && $mult == 1){
			print"*******\n$dist $evdp $phase_p $phase_s\n";
			print SAC "rh $_\n";
			print SAC "ch t1 $tp\nch t2 $ts\nch user0 $mag\n"; #Now I read location and magnitude of the reference from sac header.
			print SAC "wh\n";
			my $station =  "$kstnm.$kcmpnm";
			printf FL "%10s %10.3f %10.4f/%.4e %d %s\n",$station,$tp,$dt_dgc_p,$dt_dh_p,1,'P';
			printf FL "%10s %10.3f %10.4f/%.4e %d %s\n",$station,$ts,$dt_dgc_s,$dt_dh_s,2,'S';
		}else{
			print"*******\n$dist $evdp $phase_s\n";
			print SAC "rh $_\n";
			print SAC "ch t2 $ts\nch user0 $mag\n"; #Now I read location and magnitude of the reference from sac header.
			print SAC "wh\n";
			my $station =  "$kstnm.$kcmpnm";
			printf FL "%10s %10.3f %10.4f/%.4e %d %s\n",$station,$ts,$dt_dgc_s,$dt_dh_s,2,'S';
		}
	}
	
}
	close(FL);
}
print SAC "q\n";
close(SAC);
unlink $jk;

sub taup{
	my($model,$evdp,$dist,$phase)=@_;
	my $temp = "temp";
	system("taup_time -mod $model -h $evdp -km $dist -ph $phase > $temp"); #See TauP.
	open(TP,"<$temp");
	my @par = <TP>;
	close(TP);
	my $t0 = -100;
	shift(@par);shift(@par);shift(@par);shift(@par);shift(@par);
	my ($jk1,$jk2,$jk3,$p,$takeoff,$incident,$jk4,$jk5,$jk6);
    ($jk1,$jk2,$jk3,$t0,$p,$takeoff,$incident,$jk4,$jk5,$jk6) = split(" ",shift(@par));
	my $dtdgc = $p; #Ray parameter, horizontal slowness
	$takeoff = $takeoff*3.1415926/180;
    $t0 = sprintf("%.2f",$t0); #We have to change it to number of points in the code and make sure it can be devided by your sampling interval exactly.
    my $dtdh = - ($p/111.11)/(sin($takeoff)/cos($takeoff)); #depth slowness
	my @slow = ($t0,$dtdgc,$dtdh);
	unlink $temp;
	return (@slow);
}
