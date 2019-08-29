/*	Calculate time shift between the detected location and the reference location
 *	Miao Zhang	02/07 2015	USTC(zhmiao@mail.ustc.edu.cn)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define D2R  .017453292519943295769237
#define R2D    57.2957795130823208768
int main(int argc, char **argv) {
    int i,error;
    double lat,lon,h,Dt_Dgc,Dt_Dh,Dt;
    double x1,yp1,z1, x2,y2,z2, GCarc;
    double evla0,evlo0,evdp0,GCarc0,stla,stlo;
    error = 0;
    for(i=1; !error && i<argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
            case 'D':
                sscanf(&argv[i][2],"%lf/%lf",&Dt_Dgc,&Dt_Dh);
                break;
            case 'L':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&evla0,&evlo0,&evdp0);
                break;
            case 'E':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&lat,&lon,&h);
                break;
            case 'S':
                sscanf(&argv[i][2],"%lf/%lf",&stla,&stlo);
                break;
            default:
                error = 1;
                break;
            }
        }

    }
    if(argc < 3 || error == 1) {
        fprintf(stderr,"Usage: SHIFT -L(evla0/evlo0/evdp0) -E(lat/lon/evdp) -D(DT_Dgc/DT_Dh) -S(stla/stlo)\n");
        fprintf(stderr,"-L: reference location.\n");
        fprintf(stderr,"-E: detected location.\n");
        fprintf(stderr,"-D: horizontal slowness and vertical slowness.\n");
        fprintf(stderr,"-S: station location.\n");
        return -1;
    }

    GCinit(stla,stlo,lat,lon,&x1,&yp1,&z1,&x2,&y2,&z2,&GCarc);
    GCinit(stla,stlo,evla0,evlo0,&x1,&yp1,&z1,&x2,&y2,&z2,&GCarc0);
    Dt = Dt_Dgc*(GCarc-GCarc0) + Dt_Dh*(h-evdp0);
    printf("Dt= %lf\n",Dt);
    return 0;
}


int GCinit(lat1,lon1,lat2,lon2,x1,yp1,z1,x2,y2,z2,GCarc)
double lat1,lat2,               /*IN: (lat,lon) in degrees */
       lon1,lon2,
       *x1,*yp1,*z1,                 /* OUT: xyz coordinates */
       *x2,*y2,*z2,                  /* and  distance in radius */
       *GCarc;
{
    double
    the1,phe1,
         the2,phe2;
    the1=(90.0-lat1)*D2R;         /* convert to radius */
    phe1=lon1*D2R;
    the2=(90.0-lat2)*D2R;
    phe2=lon2*D2R;

    *x1=sin(the1)*cos(phe1);
    *yp1=sin(the1)*sin(phe1);
    *z1=cos(the1);
    *x2=sin(the2)*cos(phe2);
    *y2=sin(the2)*sin(phe2);
    *z2=cos(the2);
    *GCarc=acos((*x1)*(*x2)+(*yp1)*(*y2)+(*z1)*(*z2));
    if( fabs(*GCarc-M_PI) <=1.e-16) {
        fprintf(stderr," The great circle is not determined!\n");
        return(-1);
    }
    /*if(*GCarc<=1.e-16)  {
      fprintf(stderr,"Two same points. Program exits!\n");
      return(-1);
    }
    */
    *GCarc *= R2D;
    return(1);
}


