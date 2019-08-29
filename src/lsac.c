/********************************************************
*	Show SAC data point at specified time
*	Usage:
*		lsac time sac_files ...
*		(by Liupei Zhu)
********************************************************/

#include <stdio.h>
#include "sac.h"

int main(int argc, char **argv) {
  SACHEAD	hd;
  int		i;
  long		n;
  float		*ar;
  float		t;

  if (argc < 3) {
     fprintf(stderr, "Usage: lsac time sac_files ...\n");
     return 1;
  }

  sscanf(argv[1],"%f",&t);
  for (i=2;i<argc;i++) {
     if ((ar = read_sac(argv[i],&hd)) == NULL) {
	return -1;
     }
     n= (int) ((t-hd.b)/hd.delta);
     if (n < 0 || n > hd.npts-2)
        fprintf(stderr, "%s time out of range\n",argv[i]);
     else
	printf("%s %f\n", argv[i],ar[n]+(ar[n+1]-ar[n])*(t-hd.b-n*hd.delta)/hd.delta);
  }

  return 0;
}
