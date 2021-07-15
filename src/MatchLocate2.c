/******************************************************************************************************
 *      Match&Locate.c: Detect and locate small events from continuous seismic waveforms using templates.
 *      Copyright (c) 2013-2020 by M. Zhang
 *
 *      Method description (see Zhang and Wen, 2015a):
 *      Compared with current methods of small event detection (template matching/matched filter), 
 *      the M&L method places event detection to a lower magnitude level and extends the capability of
 *      detecting small events that have large distance separations from the template. The method has 
 *      little dependence on the accuracy of the velocity models used, and, at the same time, provides 
 *      high-precision location information of the detected small-magnitude events.
 *
 *      References:
 *      Zhang and Wen, An effective method for small event detection: Match and Locate (M&L), GJI, 2015a
 *      Zhang and Wen, Seismological Evidence for a Low‐Yield Nuclear Test on 12 May 2010 in North Korea, SRL, 2015b
 *      Zhang and Wen, Earthquake characteristics before eruptions of Japan's Ontake volcano in 2007 and 2014, GRL, 2015c
 *
 *      Usage:
 *          see the usage below.
 *
 *      Author: Miao Zhang, Dalhousie University (since April, 2019), miao.zhang@dal.ca
 *
 *      Revison History
 *          11/18 2013  M. Zhang    Initial coding
 *          06/19 2014  M. Zhang    Paper submitted
 *          02/07 2015  M. Zhang    Version 1.0 released (USTC: http://222.195.83.195/wen/codes/matchlocate)
 *          08/29 2019  M. Zhang    Version 2.0 released (Github: https://github.com/Dal-mzhang/MatchLocate2)
 *          03/13 2020  M. Zhang    Version 2.1 add bandpass filter
 *******************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <libiomp/omp.h>
#include <unistd.h>
#include "sac.h"
#define M_PI 3.14159265358979323846264338327950288
#define MAX(A,B) ((A)>(B))?(A):(B)
#define MIN(A,B) ((A)>(B))?(A):(B)
#define D2R  .017453292519943295769237
#define R2D    57.2957795130823208768
#define OUTFILE "EventCase.out"
#define CCOUT "CCdistribution.out"
#define ERROR_MAX  256
#define SAC_BANDPASS          "BP"
#define SAC_BUTTERWORTH       "BU"

#define NPEVE 9000
//max possible determined events for all locaiton.
//(!!!! Now I replace it by calculating the searching area and the time length of the continuous seismograms.)

//Define sturcture to store potential events.
typedef struct event {
    double eventLat;
    double eventLon;
    double eventH;
    double eventTime;
    double eventCoef;
    double eventNMAD;
} EVENT;


//Global parameters
double INTD;
double window,before,after;
int go,gnmax,ntrace;
float NMAD,MAD,median,THRESH;
double *Dt_Dgc,*Dt_Dh,*tshift;
int *MARKT;
double **coefsum;
double *evla0,*evlo0,*evdp0,*mag0,*gstla,*gstlo;
char **traces,**template;
double tb;
double delta=0.01,delay;

//Subroutines
void Checkdata(int, char**, int*, float*);
void Subcc(float*,float*,int,double*);
void GCinit(const double*,const double*,const double*,const double*,double*);
void Threshold(SACHEAD,double,double,double,float*,EVENT*,long int*,float);
void DetermineEvent(EVENT*,long int);
void Shiftstack(double,double,double,EVENT*,long int*);
int Time_compare_Time(const void *a,const void *b);
int Time_compare_Coef(const void *a,const void *b);
int Time_compare_NMAD(const void *a,const void *b);
int Time_compare_F(const void *a,const void *b);
int Time_compare_D(const void *a,const void *b);
int compare(const void *a,const void *b);
void Outputlast(int,double,double,double,float,float);
float Mag(float,double,double,double);
float CalculateMedian(float *,int);
float Threshold_detection(float *, int);
void taper(float *,int,float,float);
void rtrend(float *, int);
void bpcc(float *, SACHEAD, float, float);
void xapiir(float*,int,char *,double,double,int,char *,double,double,double,int);

//Main function
int main(int argc, char **argv) {
    int i,j,k,npp,error,tmark;
    SACHEAD hd,hd0,hd1;
    FILE *fp1,*fp2;
    char inputfile[256],phase[8];
    double maxlat,maxlon,maxh;
    double refevla,refevlo,refevdp;
    double lat0,lon0,h0,dlat,dlon,dh,lat,lon,h;
    long int *gk;
    extern double INTD;
    extern double window,before,after;
    extern double tb, delay,delta;
    extern int *MARKT;
    EVENT **loc1;
    EVENT *loc;
    long int *nn;
    float tleng,ppp;
    long int pn,m;
    int nlat,nlon,ndep;
    float **ref,**ref1,**obj;
    float low, high;

    error = 0;
    ntrace = 0;

    for(i=1; !error && i<argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
            case 'R':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&maxlat,&maxlon,&maxh);
                break;
            case 'I':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&dlat,&dlon,&dh);
                break;
            case 'F':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&refevla,&refevlo,&refevdp);
                break;
            case 'H':
                sscanf(&argv[i][2],"%f/%f",&THRESH,&NMAD);
                break;
            case 'T':
                sscanf(&argv[i][2],"%lf/%lf/%lf",&window,&before,&after);
                break;
            case 'D':
                sscanf(&argv[i][2],"%lf",&INTD);
                break;
            case 'B':
                sscanf(&argv[i][2],"%f/%f",&low,&high);
                break;
            case 'O':
                sscanf(&argv[i][2],"%d",&go);
                //case 0, don't output any stacked cross-correlograms.
                //case 1, output the stacked cross-correlograms of the determined events. [default]
                //case 2, output the whole stacked cross-correlograms for one searching grid.
                //case 3, output the CC value above the threshold for each searching grid.
                break;
            default:
                error = 1;
                break;
            }
        }

    }

    if(argc < 10 || error == 1) {
        fprintf(stderr, "Usage: MatchLocate2 -F(refevla/refevlo/refevdp) -R(maxlat/maxlon/maxh) -I(dlat/dlon/dh)\n");
        fprintf(stderr, "       -T(window/before/after) -H(CC/N(*MAD)) -D(INTD) -B(low/high) -O(ouput) INPUT.in\n");
        fprintf(stderr, "-F: searching center (e.g., 37.799/139.998/7.8).\n");
        fprintf(stderr, "-R: searching area (e.g., 0.05/0.05/5.0).\n");
        fprintf(stderr, "-I: searching interval (e.g., 0.01/0.01/1.0).\n");
        fprintf(stderr, "-T: time length of the reference phase (e.g., 4.0/1.0/3.0).\n");
        fprintf(stderr, "-H: cross-correlation thresholds CC && N(*MAD) (e.g., 0.3/10.0, or 0.3/0.0, or 0.0/10.0).\n");
        fprintf(stderr, "-D: keep one event within INTD sec (e.g., 6.0).\n");
        fprintf(stderr, "-B: bandpass fitering for both templates and traces (e.g., 2/8).\n");
        fprintf(stderr, "-O: output (1,2,3) or don't output (0) the cross-correlogram or CC coefficient.\n");
        fprintf(stderr, "INPUT.in: directories of templates and continuous data, horizontal and vertical slowness, etc.\n");
        return -1;
    }

    //Read the input file
    strcpy(inputfile,argv[9]);
    if((fp1 = fopen(inputfile,"r")) == NULL) {
        fprintf(stderr,"Unable to open file [INPUT.in] %s\n",inputfile);
        exit(-1);
    }

    fscanf(fp1,"%d",&ntrace);
    if(ntrace <= 0) {
        fprintf(stderr,"Please check your %s!! ntrace = %d <= zero!!\n",inputfile,ntrace);
        exit(-1);
    }

    template = (char **)malloc(sizeof(char*)*ntrace);
    traces = (char **)malloc(sizeof(char*)*ntrace);
    for(i=0; i<ntrace; i++) {
        template[i]=(char *)malloc(sizeof(char)*256);
        traces[i]=(char *)malloc(sizeof(char)*256);
    }

    Dt_Dgc =  (double *)malloc(sizeof(double)*ntrace);
    Dt_Dh  =  (double *)malloc(sizeof(double)*ntrace);
    MARKT  =  (int *)malloc(sizeof(int)*ntrace);
    tshift =  (double *)malloc(sizeof(double)*ntrace);
    gstla  =  (double *)malloc(sizeof(double)*ntrace);
    gstlo  =  (double *)malloc(sizeof(double)*ntrace);
    nn     =  (long int *)malloc(sizeof(long int)*ntrace);
    
    evla0  =  (double *)malloc(sizeof(double)*ntrace);
    evlo0  =  (double *)malloc(sizeof(double)*ntrace);
    evdp0  =  (double *)malloc(sizeof(double)*ntrace);
    mag0  =  (double *)malloc(sizeof(double)*ntrace);

    for(i=0; i<ntrace; i++) {
        fscanf(fp1, "%s %s %lf/%lf %d %s",template[i],traces[i],&Dt_Dgc[i],&Dt_Dh[i],&MARKT[i],phase);
    }
    fclose(fp1);

    //Check the data
    Checkdata(ntrace, traces, &npp, &tleng);
    //!!!!Here delay =  delta.
    delay = delta;

    //Read ref and obj
    coefsum = (double **)calloc(ntrace,sizeof(double *));
    for(i=0; i<ntrace; i++)coefsum[i] = (double *)calloc(npp,sizeof(double));
    ref = (float **)calloc(ntrace,sizeof(float *));
    for(i=0; i<ntrace; i++)ref[i] = (float *)calloc(npp,sizeof(float));
    ref1 = (float **)calloc(ntrace,sizeof(float *));
    for(i=0; i<ntrace; i++)ref1[i] = (float *)calloc(npp,sizeof(float));
    obj = (float **)calloc(ntrace,sizeof(float *));
    for(i=0; i<ntrace; i++)obj[i] = (float *)calloc(npp,sizeof(float));

    for(i=0; i<ntrace; i++) {
        	tmark=MARKT[i];
	if((ref[i]=read_sac2(template[i],&hd0,tmark,-before,after))==NULL) {
            fprintf(stderr,"Can't open template file %s\n",template[i]);
            exit(-1);
        }

        //tshift[i] = hd0.t1;
        //Mark phase arrival as tmark in SAC header 
        //Make sure it can be devided by sampling interval exactly (e.g., precision "%.1f","%.2f")
        tshift[i] = *((float *)&hd0+10+tmark);
        gstla[i] = hd0.stla;
        gstlo[i] = hd0.stlo;
	
	    evla0[i] = hd0.evla;
	    evlo0[i] = hd0.evlo;
	    evdp0[i] = hd0.evdp;
	    mag0[i] = hd0.user0; // It should be hd0.mag. There is no mag parameter in sacio.c.
        
	if((obj[i] = read_sac(traces[i],&hd))==NULL) {
            fprintf(stderr,"Can't open trace file %s\n",traces[i]);
            exit(-1);
        }
        
        //if(strcmp(hd0.kstnm,hd.kstnm)) {
            //fprintf(stderr,"The template and trace are not the same station or component %s %s\n",hd0.kstnm,hd.kstnm);
            //exit(-1);
        //}
        if(hd0.delta != hd.delta) {
            fprintf(stderr,"The template and trace are not the sampling interval %.3f %.3f\n",hd0.delta,hd.delta);
            exit(-1);
        }

        nn[i] = (long int)(((hd.e - hd.b)-window)/delay-1);
        tb = hd.b;
    }

    gnmax = nn[0];
    for(i=0; i<ntrace; i++) {
        //if(nn[i]>gnmax){gnmax=nn[i];}//use the longest waveform.
        if(nn[i]<gnmax) {
            gnmax=nn[i];   //use the shortest waveform.
        }
    }
    
    if(low > 0 && high > 0) {
        for(i=0; i<ntrace; i++) {
            //bandpass filter, eliminate the taper effect by including 0.5 sec more data
            ref1[i]=read_sac2(template[i],&hd1,tmark,-before-0.5,after+0.5);
            bpcc(ref1[i], hd1, low, high);
            bpcc(obj[i], hd, low, high);
            for(k=0;k<hd0.npts;k++){
                ref[i][k] = ref1[i][k+(int)(0.5/hd0.delta+0.5)];
            }
        }
    }

    //Compute CC
    fprintf(stderr,"********\nStart Scc!\n");
    #pragma omp parallel for firstprivate(ntrace,nn,obj,ref) private(i,k)
    for(k=0; k<gnmax; k++) {
        for(i=0; i<ntrace; i++) {
            Subcc(obj[i], ref[i], k, &coefsum[i][k]);
        }
    }
    #pragma omp barrier
    fprintf(stderr,"End Scc!\n********\n");

    for(i=0; i<ntrace; i++) {
        free(ref[i]);
        free(ref1[i]);
        free(obj[i]);
    }
    free(ref);
    free(ref1);
    free(obj);


    nlat = (int)(2*maxlat/dlat + 1);
    nlon = (int)(2*maxlon/dlon + 1);
    ndep = (int)(2*maxh/dh + 1);

    //#option 1
    //manually operation in memory setting.
    //pn = NPEVE;

    //#option 2
    //automatic memory setting.
    //!!!!Default: average two events occur in one second.
    //If it's not enough, please increase your thresholds,
    //or decrease your sampling interval, data length.
    //pn = (int)(tleng*2);
    pn = (int)(tleng*2*0.01/delta); //default sampling rate is 100 Hz.
    gk = (long int *)malloc(sizeof(long int)*nlat*nlon*ndep);
    loc = (EVENT *)malloc(sizeof(EVENT)*nlat*nlon*ndep*pn);

    loc1 = (EVENT **)malloc(sizeof(EVENT)*nlat*nlon*ndep);
    for(i=0; i<nlat*nlon*ndep; i++) {
        loc1[i] = (EVENT *)malloc(sizeof(EVENT)*pn);
    }

    ppp = nlat*nlon*ndep*pn*sizeof(EVENT)/(1024*1024);
    fprintf(stderr,"Memory: %.2lf MB.\n",2*ppp);
    if(loc == NULL || loc1 == NULL) {
        fprintf(stderr,"Memory is not enough for stucture event *loc!\nPlease shorten your data or reduce your searching area/grid.\n");
        exit(-1);
    }

    lat0 = refevla - maxlat;
    lon0 = refevlo - maxlon;
    h0 =  refevdp - maxh;

    fprintf(stderr,"Start shift and stack!\n");
    //Shift CCs and stacked them based on each potential locaiton.
    #pragma omp parallel for shared(loc1,gk) firstprivate(nlat,nlon,ndep,dlat,dlon,dh) private(lat,lon,h,m,i,j,k)
    for(m=0; m<nlat*nlon*ndep; m++) {
        i = (int)(m/(nlon*ndep));
        j = (int)((m - i*nlon*ndep)/ndep);
        k = m - i*nlon*ndep - j*ndep;
        lat = lat0 + i*dlat;
        lon = lon0 + j*dlon;
        h = h0 + k*dh;
        Shiftstack(lat,lon,h,loc1[m],&gk[m]);
    }
    #pragma omp barrier
    fprintf(stderr,"End shift and stack!\n********\nStart select events!\n");
    
    //Merge all potential events into 'loc'.
    m = 0;
    for(i=0; i<nlat*nlon*ndep; i++) {
        for(j=0; j<gk[i]; j++) {
            loc[m].eventTime = loc1[i][j].eventTime;
            loc[m].eventCoef = loc1[i][j].eventCoef;
            loc[m].eventNMAD = loc1[i][j].eventNMAD;
            loc[m].eventLat  = loc1[i][j].eventLat;
            loc[m].eventLon  = loc1[i][j].eventLon;
            loc[m].eventH    = loc1[i][j].eventH;
            m++;
        }
    }
 
    //Output the max CC at each position in order to plot the CC distribution
    //set the potential origin time to be zero.
    //12/29/2016
    if(go == 3){
	    fp2 = fopen(CCOUT,"w");
    	    if(fp1 == NULL) {
        	    fprintf(stderr,"Can't open the CC output file\n");
        	    exit(-1);
   	    }
	    fprintf(fp2,"#Time Lat. Lon. Depth Coef. N(*MAD)\n");
    	
	    for(i=0; i<nlat*nlon*ndep; i++) {
		    for(j=0; j<gk[i]; j++) {
		    fprintf(fp2,"%.3f %.4f %.4f %.4f %.6f %.4f\n",loc1[i][j].eventTime,loc1[i][j].eventLat,loc1[i][j].eventLon,loc1[i][j].eventH,loc1[i][j].eventCoef,loc1[i][j].eventNMAD);
		    }

        }
	    fclose(fp2);
    }

    //Determine events and output the stacked cross-correlograms [go = 1].
    DetermineEvent(loc,m);
    fprintf(stderr,"End select events!\n");

    //Free memory
    free(loc);
    free(gk);
    for(i=0; i<nlat*nlon*ndep; i++) free(loc1[i]);
    free(loc1);
    for(i=0; i<ntrace; i++) {
        free(template[i]);
        free(traces[i]);
        free(coefsum[i]);
    }
    free(template);
    free(traces);
    free(coefsum);
    free(Dt_Dgc);
    free(Dt_Dh);
    free(MARKT);
    free(tshift);
    free(gstla);
    free(gstlo);
    free(nn);
    free(evla0);
    free(evlo0);
    free(evdp0);
    free(mag0);
    return 0;
}

//Here just compute the normalized CC value. (time lag at t = 0)
void Subcc(float *obj, float *ref,  int k, double *fp) {
    int j;
    extern double tb;
    extern double delta,delay,window;
    double cc,norm,normMaster,cczero;
    double tstart2,tend2;
    int n1,n2;

    tstart2 = tb + k*delay;
    n1 = (int)((tstart2-tb)/delta+0.5);
    tend2 = tstart2 + window;
    n2 = (int)((tend2-tb)/delta+0.5);
    cc = 0.0;
    norm = 0.0;
    normMaster = 0.0;
    for(j=0; j<(n2-n1); j++) {
        cc += ref[j]*obj[n1+j];
        norm += obj[n1+j]*obj[n1+j];
        normMaster += ref[j]*ref[j];
    }
    if(norm == 0 || normMaster == 0){*fp = 0.0; return;}
    cczero = cc/(sqrt(norm)*sqrt(normMaster));
    //cczero = cc/(sqrt(norm+1.0e-20)*sqrt(normMaster+1.0e-20));
    *fp = cczero;
}

//shift cross-correlograms and stacked them based on the potential locations.
void Shiftstack(double lat, double lon, double h, EVENT *loc2,long int *KK) {
    extern int ntrace,gnmax,go;
    extern double *Dt_Dgc,*Dt_Dh,*tshift;
    extern double *evla0,*evlo0,*evdp0,*gstla,*gstlo,**coefsum;
    extern double delay,before;
    double **coefout;
    float *sumout;
    float tshiftmax,threshold0;
    double Dt,GCarc,GCarc0;
    SACHEAD hd1;
    int i,k,nshift,gnccmax;
    char outfile[128];

    sumout = (float *)calloc(gnmax,sizeof(float));
    coefout = (double **)calloc(ntrace,sizeof(double *));
    for(i=0; i<ntrace; i++) {
        coefout[i] = (double *)calloc(gnmax,sizeof(double));
    }
    
    tshiftmax = 0.0;
    for(i=0; i<ntrace; i++) {
        GCinit(&gstla[i],&gstlo[i],&lat,&lon,&GCarc);
        GCinit(&gstla[i],&gstlo[i],&evla0[i],&evlo0[i],&GCarc0);
        Dt = Dt_Dgc[i]*(GCarc-GCarc0) + Dt_Dh[i]*(h-evdp0[i]);
        nshift = (int)((tshift[i] + Dt)/delay+0.5);

        for(k=0; k<gnmax; k++) {
            if((k+nshift)<gnmax && (k + nshift)>=0) {
                //Here ignore several seconds (tshift)seismograms at the beginning.
                coefout[i][k] = coefsum[i][k+nshift];
            }
            else {
                coefout[i][k]=0.0;
            }
        }
        if(tshift[i] > tshiftmax){tshiftmax = tshift[i];}
    }

    //Remove the blank Cross-correlogram (segment of travel time t1)
    gnccmax = gnmax - (int)(tshiftmax/delay);   
    for(k=0; k<gnmax; k++) {
        for(i=0; i<ntrace; i++) {
            sumout[k] += coefout[i][k];
        }
        sumout[k] /= ntrace;
    }


    hd1 = sachdr(delay,gnccmax,tb+before);
    hd1.evla = lat;
    hd1.evlo = lon;
    hd1.evdp = h;

    //Output stacked cross-correlograms for each searching grid.
    if(go == 2 ){
        sprintf(outfile,"%07.4f_%08.4f_%06.2f.stack",lat,lon,h);
        write_sac(outfile,hd1,sumout);
    }
    

    threshold0 = Threshold_detection(sumout, hd1.npts);

    //Determine potential events.
    Threshold(hd1,lat,lon,h,sumout,loc2,KK,threshold0);

    free(sumout);
    for(i=0; i<ntrace; i++) free(coefout[i]);
    free(coefout);
}

//The coef. threshold
void Threshold(SACHEAD hd,double evla, double evlo, double evdp, float *ar, EVENT *loc2, long int *KK, float threshold) {
    long int j,k;
    extern double delay;
    extern float median, MAD, THRESH;

    k = 0; 
    for(j=0; j<hd.npts; j++) {
        if(ar[j] >= threshold && ar[j] >= THRESH) {
            loc2[k].eventTime = (double)(hd.b) + delay*j;
            loc2[k].eventCoef = ar[j];
            loc2[k].eventNMAD = (ar[j] - median)/MAD;
            loc2[k].eventLat = evla;
            loc2[k].eventLon = evlo;
            loc2[k].eventH = evdp;
            //Here is a test
            //fprintf(stderr,"%ld %lf %f %f %f %f %f\n",k,loc2[k].eventTime,loc2[k].eventCoef,loc2[k].eventNMAD,hd.evla,hd.evlo,hd.evdp);
            k++;
        }
    }
    *KK = k;
}

//Select events
void DetermineEvent(EVENT *loc,long int NN) {
    EVENT *loc2;
    int i,j,k,nn;
    FILE *fp;
    float mag;
    extern float NMAD;

    loc2 = (EVENT *)malloc(sizeof(EVENT)*NN);
    if(loc2 == NULL) {
        fprintf(stderr,"Can't locate loc2 in DetermineEvent\n");
        exit(-1);
    }
    
    //Sort by coefficient or MAD (from large to small).
    if(NMAD > 1){
        qsort(loc,NN,sizeof(EVENT),Time_compare_NMAD);
    }else{
        qsort(loc,NN,sizeof(EVENT),Time_compare_Coef);
    }

    fp = fopen(OUTFILE,"w");
    if(fp == NULL) {
        fprintf(stderr,"Can't open output file in DetermineEvent\n");
        exit(-1);
    }
    fprintf(fp,"#Event   Time      Lat.      Lon.        Depth    Mag.    Coef.      N(*MAD)  NumTrace\n");

    if(NN > 0) {
        loc2[0].eventTime = loc[0].eventTime;
        loc2[0].eventCoef = loc[0].eventCoef;
        loc2[0].eventNMAD = loc[0].eventNMAD;
        loc2[0].eventLat = loc[0].eventLat;
        loc2[0].eventLon = loc[0].eventLon;
        loc2[0].eventH = loc[0].eventH;

        k = 1;
        for(i=1; i<NN; i++) {
            nn = 0;
            for(j=0; j<k; j++) {
                if(fabs(loc[i].eventTime - loc2[j].eventTime) > INTD) {
                    nn++;
                }
            }

            if(nn == k || loc[i].eventCoef == 1.00) {
                loc2[k].eventTime = loc[i].eventTime;
                loc2[k].eventCoef = loc[i].eventCoef;
                loc2[k].eventNMAD = loc[i].eventNMAD;
                loc2[k].eventLat = loc[i].eventLat;
                loc2[k].eventLon = loc[i].eventLon;
                loc2[k].eventH = loc[i].eventH;
                k++;
            }

        }

        //Sort by origin time (from small to large).
        qsort(loc2,k,sizeof(EVENT),Time_compare_Time);

        for(i=0; i<k; i++) {
            mag = Mag(loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH);
            fprintf(fp,"%4d   %9.3lf   %7.4f   %8.4f   %6.2f   %5.2f    %6.4f    %7.4f   %4d\n",i+1,
                    loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,mag,loc2[i].eventCoef,loc2[i].eventNMAD,ntrace);
            if(go == 1) {
                Outputlast(i+1,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,loc2[i].eventTime,loc2[i].eventCoef);
            }
        }
    }
    fclose(fp);
    free(loc2);
}

//calculate median value
float CalculateMedian(float *arrValue, int max)
{
    float median = 0;
    float *value;
    int i, j;
    float temp;
    value = (float *)malloc(max*sizeof(float));
    for(i = 0; i < max; i++)
                value[i] = arrValue[i];
    for(i = 0; i < max; i++)
    {
        for(j = 0; j < max - i - 1; j++)
        {
            if(value[j] > value[j + 1])
            {
                temp = value[j];
                value[j] = value[j + 1];
                value[j + 1] = temp;
            }
        }
    }
    if( (max % 2) == 1)
    {
        median =  value[ (max + 1) / 2 - 1];
    }
    else
    {
        median = (value[max / 2] + value[max / 2 - 1]) / 2;
    }
    free(value);
    return median;
}

//calculate the threshold of detection
float Threshold_detection(float *cc_sum, int npts){
    extern float median, MAD, NMAD;
    float *cc_sum2,*cc_sum1;

    cc_sum1=(float *)calloc(npts,sizeof(float ));
    for(int i=0;i<npts;i++)cc_sum1[i]=cc_sum[i];
    cc_sum2=(float *)calloc(npts,sizeof(float ));
    qsort(cc_sum1, npts, sizeof(float), compare);
    if(npts%2 == 0)median = (cc_sum1[npts/2]+cc_sum1[npts/2-1])/2;
    if(npts%2 == 1)median = cc_sum1[(npts-1)/2];
    for(int i=0; i<npts; i++)cc_sum2[i] = fabsf(cc_sum1[i]-median);
    qsort(cc_sum2, npts, sizeof(float), compare);
    if(npts%2 == 0)MAD = (cc_sum2[npts/2]+cc_sum2[npts/2-1])/2;
    if(npts%2 == 1)MAD = cc_sum2[(npts-1)/2];
    free(cc_sum1);
    free(cc_sum2);
    return median+NMAD*MAD;
}

//Determine the magnitude
float Mag(float initime,double lat, double lon, double h) {
    double t1,t2,Dt;
    float *master,*ar,*ratio,*mag;
    double GCarc,GCarc0;
    extern int ntrace,*MARKT;
    extern double *gstla,*gstlo,*Dt_Dgc,*Dt_Dh,*mag0,before,window,after;
    int i,k,npts1,npts2,tmark;
    float armax,mastermax,median;
    SACHEAD hd0,hd1;
    int tmark2 = -3;

    ratio = (float *)malloc(ntrace*sizeof(float));
    mag = (float *)malloc(ntrace*sizeof(float));
    for(i=0; i<ntrace; i++) {
	    tmark=MARKT[i];
        if( (master=read_sac2(template[i],&hd0,tmark,-before,after)) == NULL) {
            fprintf(stderr,"erro master in Mag\n");
            exit(-1);
        }
        npts1 = hd0.npts;
        mastermax = master[0];
        for(k=0; k<npts1; k++) {
            if(fabs(master[k]) > mastermax)mastermax = fabs(master[k]);
        }

        GCinit(&gstla[i],&gstlo[i],&lat,&lon,&GCarc);
        GCinit(&gstla[i],&gstlo[i],&evla0[i],&evlo0[i],&GCarc0);
        Dt = Dt_Dgc[i]*(GCarc-GCarc0) + Dt_Dh[i]*(h-evdp0[i]);
        t1 = initime + tshift[i] + Dt - before;
        t2 = t1 + window;
        if( (ar=read_sac2(traces[i],&hd1,tmark2,t1,t2)) == NULL) {
            fprintf(stderr,"erro event in Mag\n");
            exit(-1);
        }
        npts2 = hd1.npts;
        armax = ar[0];
        for(k=0; k<npts2; k++) {
            if(fabs(ar[k]) > armax)armax = fabs(ar[k]);
        }
	if(armax == 0 || mastermax == 0) continue;
        ratio[i] = armax/mastermax;
    	mag[i] = mag0[i] + log10(ratio[i]);
    }
    free(master);
    free(ar);
    median = CalculateMedian(mag,ntrace);
    free(ratio);free(mag);
    return median;
}

//Based on the determined location, output the stacked ccs.
void Outputlast(int ll,double lat, double lon, double h,float initime,float coef) {
    extern int ntrace,gnmax,go;
    extern double *Dt_Dgc,*Dt_Dh,*tshift,**coefsum,tb;
    extern double *evla0,*evlo0,*evdp0,*gstla,*gstlo;
    extern double delay,INTD,before;
    double **coefout;
    float *sumout;
    double Dt,GCarc,GCarc0;
    SACHEAD hd1;
    int i,k,nshift;
    char temp[8],outfile[128];

    sumout = (float *)calloc(gnmax,sizeof(float));
    coefout = (double **)calloc(ntrace,sizeof(double *));
    for(i=0; i<ntrace; i++) {
        coefout[i] = (double *)calloc(gnmax,sizeof(double));
    }
    for(i=0; i<ntrace; i++) {
        GCinit(&gstla[i],&gstlo[i],&lat,&lon,&GCarc);
        GCinit(&gstla[i],&gstlo[i],&evla0[i],&evlo0[i],&GCarc0);
        Dt = Dt_Dgc[i]*(GCarc-GCarc0) + Dt_Dh[i]*(h-evdp0[i]);
        nshift = (int)((tshift[i] + Dt)/delay+0.5);
        for(k=0; k<gnmax; k++) {
            if((k+nshift)<gnmax && (k + nshift)>=0) {
                coefout[i][k] = coefsum[i][k+nshift];
            }
            else {
                coefout[i][k]=0.0;
            }
        }
    }
    for(k=0; k<gnmax; k++) {
        for(i=0; i<ntrace; i++) {
            sumout[k] += coefout[i][k];
        }
        sumout[k] /= ntrace;
    }

    /*
    // output all time stack
        hd1 = sachdr(delay,gnmax,tb+before);
        hd1.evla = lat;
        hd1.evlo = lon;
        hd1.evdp = h;
        hd1.t1 = initime;
        sprintf(temp,"%4.2f",coef);
        strcpy(hd1.kt1,temp);
        sprintf(outfile,"%04d_%08.4f_%07.4f_%08.4f_%06.2f.stack",ll,initime,lat,lon,h);
        write_sac(outfile,hd1,sumout);
    */
//  output just [t-INTD, t1+INTD];
    double timeb = initime - INTD;
    double timee = initime + INTD;
    int n1, n2;
    float *sumcut;
    n1 = (int)((timeb-tb-before)/delay+0.5);
    n2 = (int)((timee-tb-before)/delay+0.5);
    sumcut = (float *)calloc((n2-n1),sizeof(float));
    for(k=0; k<(n2-n1); k++) {
        sumcut[k] = sumout[k+n1];
    }
    hd1 = sachdr(delay,n2-n1,timeb);
    hd1.evla = lat;
    hd1.evlo = lon;
    hd1.evdp = h;
    hd1.t1 = initime;
    sprintf(temp,"%6.4f",coef);
    strcpy(hd1.kt1,temp);
    sprintf(outfile,"%04d_%08.2f_%07.4f_%08.4f_%06.2f.stack",ll,initime,lat,lon,h);
    write_sac(outfile,hd1,sumcut);
    free(sumcut);
    free(sumout);
    for(i=0; i<ntrace; i++) free(coefout[i]);
    free(coefout);
}

//Compare function
int Time_compare_Time(const void *a,const void *b) {
    double v1 = ((EVENT *)a)->eventTime;
    double v2 = ((EVENT *)b)->eventTime;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_Coef(const void *b,const void *a) {
    double v1 = ((EVENT *)a)->eventCoef;
    double v2 = ((EVENT *)b)->eventCoef;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_NMAD(const void *b,const void *a) {
    double v1 = ((EVENT *)a)->eventNMAD;
    double v2 = ((EVENT *)b)->eventNMAD;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_F(const void *a,const void *b) {
    float v1 = *(float *)a;
    float v2 = *(float *)b;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_D(const void *a,const void *b) {
    int v1 = *(int *)a;
    int v2 = *(int *)b;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int compare (const void *a, const void *b) {
    const float *da = (const float *) a;
    const float *db = (const float *) b;
    return  (*da > *db) - (*da < *db);
}

//Compute great circle (gcarc)
void GCinit(const double *lat1,const double *lon1,const double *lat2,const double *lon2,double *GCarc)
{
    double x1,yp1,z1,x2,y2,z2;
    double the1,phe1,the2,phe2;
    the1=(90.0-*lat1)*D2R;         /* convert to radius */
    phe1=(*lon1)*D2R;
    the2=(90.0-*lat2)*D2R;
    phe2=(*lon2)*D2R;

    x1=sin(the1)*cos(phe1);
    yp1=sin(the1)*sin(phe1);
    z1=cos(the1);
    x2=sin(the2)*cos(phe2);
    y2=sin(the2)*sin(phe2);
    z2=cos(the2);
    *GCarc=acos((x1)*(x2)+(yp1)*(y2)+(z1)*(z2));
    if( fabs(*GCarc-M_PI) <=1.e-16) {
        fprintf(stderr," The great circle is not determined! (GCinit)\n");
        exit(1);
    }
    *GCarc *= R2D;
}

// check the continuous data, you can ignore below part.
void Checkdata(int ntrace, char **traces, int *npp, float *tleng) {
    float *nnb,*nne;
    double *nnd;
    int *nnl, i;
    extern double delta;
    SACHEAD hd;
    nnb     =  (float *)malloc(sizeof(float)*ntrace);
    nne     =  (float *)malloc(sizeof(float)*ntrace);
    nnd     =  (double *)malloc(sizeof(double)*ntrace);
    nnl     =  (int *)malloc(sizeof(int)*ntrace);
    for(i=0; i<ntrace; i++) {
        read_sachead(traces[i],&hd);
        nnb[i] = hd.b;
        nne[i] = hd.e;
        nnd[i] = (double)(hd.delta);
        nnl[i] = (int)(((hd.e - hd.b))/hd.delta+0.5);
        *tleng = hd.e - hd.b;
        if(hd.o != 0) {
            fprintf(stderr,"Please set origin time ZERO [ch o 0] in your continuous data.\n");
            exit(-1);
        }
    }

    qsort(nnb,ntrace,sizeof(nnb[0]),Time_compare_F); // sort from small to large
    qsort(nne,ntrace,sizeof(nne[0]),Time_compare_F);
    qsort(nnd,ntrace,sizeof(nnd[0]),Time_compare_F);
    qsort(nnl,ntrace,sizeof(nnl[0]),Time_compare_D);


    if(fabs(nnb[ntrace-1] - nnb[0]) > 0.5*nnd[0]) {
        fprintf(stderr,"Please Check Your Continuous Data!!! Not the same beginning time.\n");
        exit(-1);
    }
    if(fabs(nnd[ntrace-1] - nnd[0]) > 1.0e-5) {
        fprintf(stderr,"Please Check Your Continuous Data!!! Not the same sampling invertal.\n");
        exit(-1);
    }

    *npp = nnl[0];
    if(fabs(nne[ntrace-1] - nne[0]) > 1.0e-3) {
        fprintf(stderr,"Please Check Your Continuous Data!!! Not the same end time.\n");
        fprintf(stderr,"Shortest end time: %.2f sec; Longest end time %.2f sec. The shortest one would be used.\n",
                nne[0],nne[ntrace-1]);
    }

    delta = (double)(nnd[0]);

    free(nnb);
    free(nne);
    free(nnd);
    free(nnl);
}


//hanning taper before bandpass filter
void taper(float *yarray,int nlen,float start,float end){
    float ang,cs;
    int m1,m2,m3,m4,m5;
    int i,j,k,xi;
         
	m1 = (int)(nlen*start + 0.5);
	m2 = m1 + 1;

	ang = 3.1415926/(float)(m1);
	
	for(i=0;i<=m1;i++){
	xi=i;
	cs = (1 - cos(xi*ang))/2.0;
	yarray[i] = yarray[i]*cs;
	}
	
	m3 = (int)(nlen*end + 0.5);
	m5 = nlen - m3-1;
	m4 = m5+1;
	ang = 3.1415926/(float)(m3);
	
	for(k=m2;k<=m5;k++){
	yarray[k] = yarray[k];
	}
	for(j=m4;j<nlen;j++){
	xi = j + 1 - nlen;
	cs = (1 - cos(xi*ang))/2.0;
	yarray[j] = yarray[j]*cs;
	}
}

/* remove trend a*i + b */
void rtrend(float *y, int n) {
     int i;
     double a, b, a11, a12, a22, y1, y2;
     y1 = y2 = 0.;
     for(i=0;i<n;i++) {
       y1 += i*y[i];
       y2 += y[i];
     }
     a12 = 0.5*n*(n-1);
     a11 = a12*(2*n-1)/3.;
     a22 = n;
     b = a11*a22-a12*a12;
     a = (a22*y1-a12*y2)/b;
     b = (a11*y2-a12*y1)/b;
     for(i=0;i<n;i++) {
       y[i] = y[i] - a*i - b;
     }
}

//do bandpass filtering for tempates and traces
void bpcc(float *yarray, SACHEAD hd, float low, float high) {
    /* Local variables */
    //float low, high;
    double attenuation, transition_bandwidth;
    int nlen;
    //SACHEAD hd;
    double delta_d;
    int order;
    //float *yarray;
    int passes;
    float total,sum,mean,taperb;
    int j;

    sum = 0.0;
    //rmean
    for(j=0;j<hd.npts;j++) {sum += yarray[j];}
    for(j=0;j<hd.npts;j++) {yarray[j] = yarray[j] - sum/hd.npts;}
    delta_d = hd.delta;
    nlen = hd.npts;

    //rtrend
    rtrend(yarray,nlen);
    /*taper function*/
    //taper(yarray,hd.npts,0.05,0.05); //sac default hanning window 0.05
    taperb = 0.0001;
    if(hd.npts<20000){taperb = 0.01;}
    taper(yarray,hd.npts,taperb,taperb); //taper first 50 points
    passes = 2;
    order  = 4;
    transition_bandwidth = 0.0;
    attenuation = 0.0;
    xapiir(yarray, nlen, SAC_BUTTERWORTH, transition_bandwidth, attenuation, order, SAC_BANDPASS, low, high, delta_d, passes);
    //write_sac("test.sac",hd,yarray);
}
