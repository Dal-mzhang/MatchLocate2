#An example to plot CC distribution for Match&Locate
#By Miao Zhang, Dalhousie University, miao.zhang@dal.ca
#!/usr/bin/perl -w
$ps = "CCdistribution.ps";
$CCout = "../CCdistribution.out";
$trange = 1;
$range = 0.008;
##template location (from catalog)
$latref = "37.788";
$lonref = "140.001";
$lat0 = $latref;
$lon0 = $lonref;
##slave location (from catalog)
$latslave = "37.791";
$lonslave = "140.005";

$lat1 = $lat0 - $range*1.0;
$lat2 = $lat0 + $range*1.0;
$lon1 = $lon0 - $range*1.0;
$lon2 = $lon0 + $range*1.0;

$R1 = "$lon1/$lon2/$lat1/$lat2";
$B1 = "a0.005f0.002/a0.005f0.002";
$J1 = "X5i";

`psxy -R$R1 -J$J1 -P -K /dev/null > $ps`;
system "psbasemap -R$R1 -J$J1 -B${B1}WeSn -K -O >> $ps";

open(JK1,"<$CCout");
@dat1 = <JK1>;
close(JK1);
shift(@dat1);

$ccmax = 0;
foreach $_(@dat1){
	chomp($_);
	($otime,$lat,$lon,$dep,$cc,$jk) = split(" ",$_);
	if($cc > $ccmax){$ccmax = $cc;
	##located slave event by M&L
	$latslave_ml = $lat;
	$lonslave_ml = $lon;}
}

$cc0 = $ccmax*0.86;
$cct = $ccmax*0.65;

open(GMT,"|xyz2grd -Ggrdfile -I0.0001/0.0001 -R$R1 -V");
foreach $_(@dat1){
	chomp($_);
	($otime,$lat,$lon,$dep,$cc,$jk) = split(" ",$_);
	if($cc > $cct && $otime > -$trange && $otime < $trange){print GMT "$lon $lat $cc\n";}
}
close(GMT);



$cpt = "cpt.file";
system("makecpt -Chot -T$cct/$ccmax/0.02 -I -V > $cpt");
system("grdimage grdfile -C$cpt -R$R1 -J$J1 -K -O -V >> $ps");

#kind of way to estimate uncertainty
`grdcontour -R -J  grdfile -C$cc0 -L0/$ccmax -K -O -W2p,white >> $ps`;

open(GMT,"|psxy -R$R1 -J$J1 -K -O -Sa0.15i -Gblack >> $ps");
print GMT "$lonref $latref\n";
close(GMT);
&PSTEXT($lonref,$latref-0.0005,15,0,4,MC,"Template (Catalog)",$ps);
open(GMT,"|psxy -R$R1 -J$J1 -K -O -Sa0.15i -Ggray >> $ps");
print GMT "$lonslave $latslave\n";
close(GMT);
&PSTEXT($lonslave,$latslave-0.0005,15,0,4,MC,"Slave (Catalog)",$ps);
open(GMT,"|psxy -R$R1 -J$J1 -K -O -Sa0.15i -Gblue >> $ps");
print GMT "$lonslave_ml $latslave_ml\n";
close(GMT);
&PSTEXT($lonslave_ml,$latslave_ml-0.0005,15,0,4,MC,"Slave (M&L)",$ps);

$latm = ($lat1+$lat2)/2;
$lonm = ($lon1+$lon2)/2;
&PSTEXT($lon1-0.003,$latm,20,90,4,MC,"Latitude",$ps);
&PSTEXT($lonm,$lat1-0.0015,20,0,4,MC,"Longitude",$ps);
`psscale -D6/1.5/8/0.3h -K -O -C$cpt -B0.05/:"CC": >> $ps`;

unlink "grdfile";
unlink $cpt;
`psxy -R -J -O  /dev/null >> $ps`;
`ps2raster -A -P -E500 $ps`;
sub PSTEXT{
        my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
        open(GMT,"| pstext -R -J -K -O -N >> $ps");
                print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
        close(GMT);
}

