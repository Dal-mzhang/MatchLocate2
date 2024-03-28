/******************************************************************************************************
//Miao Zhang, Dalhousie University
//miao.zhang@dal.ca
//Select best events in each day (they may be detected by multiple templtes).
 *******************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <libiomp/omp.h>
#include <math.h>
#include <unistd.h>
#include "sac.h"

#define NPEVE 100000	//max possible determined events for all locaiton

//define sturcture to store the potential events.
typedef struct event {
    double eventLat;
    double eventLon;
    double eventH;
    double eventTime;
    double eventCoef;
    double eventNMAD;
    double eventMag;
    int eventtr;
    char eventRef[128];
} EVENT;

//Globle parameter
float NMAD,INTD,THRESH;
char outfile[128];

//Subroutine
void DetermineEvent(EVENT *loc, long int);
int Time_compare_Time(const void *a,const void *b);
int Time_compare_NMAD(const void *a,const void *b);
int Time_compare_Coef(const void *a,const void *b);

//Main function begin here.
int main(int argc, char **argv) {
    int i,k,kk,error;
    FILE *fp1;
    char inputfile[128];
    extern float NMAD,INTD,THRESH;
    long int gk;
    EVENT *loc, *loc1;
    extern char outfile[128];
    float HHN;

    error = 0;
    for(i=1; !error && i<argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
            case 'D':
                sscanf(&argv[i][2],"%f",&INTD);
                break;
            case 'H':
                sscanf(&argv[i][2],"%f/%f",&THRESH,&NMAD);
                break;
            default:
                error = 1;
                break;
            }
        }

    }
    if(argc < 3 || error == 1) {
        fprintf(stderr,"Usage: SelectFinal -H(CC/NMAD) -D(INTD) Eventfile\n");
        fprintf(stderr,"-H: CC >= THRESH &&  N > NMAD (N*MAD) (e.g., 0.3/0.0, 0.0/10.0, or 0.3/10.0)\n");
        fprintf(stderr,"-D: keep one event within INTD sec (e.g., 6.0)\n");
        fprintf(stderr,"Eventfile: potential events detected by all templates.\n");
        return -1;
    }

    //Read the input file
    for(i=3; i<argc; i++) {
        strcpy(inputfile,argv[i]);
        sprintf(outfile,"%s.final",inputfile);
        if((fp1 = fopen(inputfile,"r")) == NULL) {
            fprintf(stderr,"Unable to open file [Event files] %s\n",inputfile);
            exit(-1);
        }

        loc = (EVENT *)malloc(sizeof(EVENT)*NPEVE);
        loc1 = (EVENT *)malloc(sizeof(EVENT)*NPEVE);
        if(loc == NULL || loc1 == NULL) {
            fprintf(stderr,"There are no enough memory for loc and loc1\n");
            exit(-1);
        }


        k = 0;
        kk = 0;
        //fprintf(stderr,"Start read all potential events!\n");
        while((fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %d %s", &loc1[k].eventTime, &loc1[k].eventLat,
                      &loc1[k].eventLon, &loc1[k].eventH, &loc1[k].eventMag, &loc1[k].eventCoef,
                      &loc1[k].eventNMAD, &loc1[k].eventtr, &loc1[k].eventRef)) != EOF) {
            if(loc1[k].eventNMAD >= NMAD && loc1[k].eventCoef >= THRESH) {
                loc[kk].eventTime = loc1[k].eventTime;
                loc[kk].eventLat = loc1[k].eventLat;
                loc[kk].eventLon = loc1[k].eventLon;
                loc[kk].eventH = loc1[k].eventH;
                loc[kk].eventMag = loc1[k].eventMag;
                loc[kk].eventCoef = loc1[k].eventCoef;
                loc[kk].eventNMAD = loc1[k].eventNMAD;
                loc[kk].eventtr = loc1[k].eventtr;
                strcpy(loc[kk].eventRef,loc1[k].eventRef);
                kk++;
            }
            k++;
        }

        gk = kk;
        //fprintf(stderr,"Start select events!\n");
        DetermineEvent(loc,gk);
        //fprintf(stderr,"End select events!\n");

        //Free memory
        free(loc1);
        free(loc);
    }
    return 0;
}

void DetermineEvent(EVENT *loc,long int gk) {
    EVENT *loc2;
    extern char outfile[128];
    extern float INTD;
    int i,j,k,nn;
    FILE *fp;

    loc2 = (EVENT *)malloc(sizeof(EVENT)*gk);
    if(loc2 == NULL) {
        fprintf(stderr,"Can't locate loc2\n");
        exit(-1);
    }

    qsort(loc,gk,sizeof(EVENT),Time_compare_NMAD);
    //qsort(loc,gk,sizeof(EVENT),Time_compare_Coef);

    fp = fopen(outfile,"w");
    if(fp == NULL) {
        fprintf(stderr,"Can't open output file in DetermineEvent\n");
        exit(-1);
    }
    fprintf(fp,"#Event     Time     Lat.      Lon.        Depth    Mag.    Coef.     N(*MAD)    Num_Trace     Reference\n");

    if(gk > 0) {
        loc2[0].eventTime = loc[0].eventTime;
        loc2[0].eventLat = loc[0].eventLat;
        loc2[0].eventLon = loc[0].eventLon;
        loc2[0].eventH = loc[0].eventH;
        loc2[0].eventMag = loc[0].eventMag;
        loc2[0].eventCoef = loc[0].eventCoef;
        loc2[0].eventNMAD = loc[0].eventNMAD;
        loc2[0].eventtr = loc[0].eventtr;
        strcpy(loc2[0].eventRef,loc[0].eventRef);

        k = 1;
        for(i=1; i<gk; i++) {
            nn = 0;
            for(j=0; j<k; j++) {
                if(fabs(loc[i].eventTime - loc2[j].eventTime) > INTD) {
                    nn++;
                }
            }

            if(nn == k || loc[i].eventCoef == 1.0 ) {
                loc2[k].eventTime = loc[i].eventTime;
                loc2[k].eventLat = loc[i].eventLat;
                loc2[k].eventLon = loc[i].eventLon;
                loc2[k].eventH = loc[i].eventH;
                loc2[k].eventMag = loc[i].eventMag;
                loc2[k].eventCoef = loc[i].eventCoef;
                loc2[k].eventNMAD = loc[i].eventNMAD;
                loc2[k].eventtr = loc[i].eventtr;
                strcpy(loc2[k].eventRef,loc[i].eventRef);
                k++;
            }

        }

        qsort(loc2,k,sizeof(EVENT),Time_compare_Time);

        for(i=0; i<k; i++) {
            fprintf(fp,"%4d   %9.3lf   %7.4f   %8.4f   %6.2f   %5.2f    %6.4f    %8.4f    %3d    %s\n",i+1,
                    loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,loc2[i].eventMag,
                    loc2[i].eventCoef,loc2[i].eventNMAD,loc2[i].eventtr,loc2[i].eventRef);
        }
    }
    fclose(fp);
    free(loc2);
}

//Compare function
int Time_compare_Time(const void *a,const void *b) {
    double v1 = ((EVENT *)a)->eventTime;
    double v2 = ((EVENT *)b)->eventTime;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

//Compare function
int Time_compare_NMAD(const void *b,const void *a) {
    double v1 = ((EVENT *)a)->eventNMAD;
    double v2 = ((EVENT *)b)->eventNMAD;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

//Compare function
int Time_compare_Coef(const void *b,const void *a) {
    double v1 = ((EVENT *)a)->eventNMAD;
    double v2 = ((EVENT *)b)->eventNMAD;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}
