
#ifndef __SACIO_H__
#define __SACIO_H__

/** 
 * Define an automatic string length from C
 *
 */
#define SAC_STRING_LENGTH   -1

#ifndef TRUE
#define TRUE  1
#endif /* TRUE */

#ifndef FALSE
#define FALSE 0
#endif /* FALSE */

#define SAC_FLOAT_UNDEFINED      -12345.0
#define SAC_REAL_UNDEFINED       SAC_FLOAT_UNDEFINED
#define SAC_INT_UNDEFINED        -12345
#define SAC_INTEGER_UNDEFINED    SAC_INT_UNDEFINED
#define SAC_NUMBER_UNDEFINED     SAC_INT_UNDEFINED
#define SAC_CHAR_UNDEFINED       "-12345  "
#define SAC_CHARACTER_UNDEFINED  SAC_CHAR_UNDEFINED

void getfhv(char      *kname, 
            float     *fvalue, 
            int       *nerr, 
            int        kname_s);

void getihv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s, 
            int        kvalue_s);

void getkhv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s, 
            int        kvalue_s);

void getlhv(char      *kname, 
            int       *lvalue, 
            int       *nerr, 
            int        kname_s);

void getnhv(char      *kname, 
            int       *nvalue, 
            int       *nerr, 
            int        kname_s);

void newhdr ();

void rsac1(char      *kname, 
           float      yarray[], 
           int       *nlen, 
           float     *beg, 
           float     *del, 
           int       *max_, 
           int       *nerr, 
           int        kname_s);

void rsac2(char      *kname, 
           float      yarray[], 
           int       *nlen, 
           float      xarray[], 
           int       *max_, 
           int       *nerr, 
           int        kname_s);

void rsach(char      *kname,
           int       *nerr,
           int        kname_s);

void setfhv(char      *kname, 
            float     *fvalue, 
            int       *nerr, 
            int        kname_s);

void setihv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s,
            int        kvalue_s);

void setkhv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s,
            int        kvalue_s);

void setlhv(char      *kname, 
            int       *lvalue, 
            int       *nerr, 
            int        kname_s);

void setnhv(char      *kname, 
            int       *nvalue, 
            int       *nerr, 
            int        kname_s);

void wsac0(char      *kname, 
           float     *xarray, 
           float     *yarray,
           int       *nerr,
           int        kname_s);

void wsac1(char      *kname, 
           float     *yarray, 
           int       *nlen, 
           float     *beg, 
           float     *del, 
           int       *nerr, 
           int        kname_s);

void wsac2(char      *kname, 
           float     *yarray, 
           int       *nlen, 
           float     *xarray, 
           int       *nerr, 
           int        kname_s);

void wsac3(char      *kname, 
           float     *xarray, 
           float     *yarray,
           int       *nerr,
           int        kname_s);

void getbbv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s, 
            int        kvalue_s);

void writebbf(char      *kname, 
              int       *nerr, 
              int        kname_s);

void readbbf(char      *kname, 
             int       *nerr, 
             int        kname_s);


void setbbv(char      *kname, 
            char      *kvalue, 
            int       *nerr, 
            int        kname_s, 
            int        kvalue_s);



#endif /* __SACIO_H__ */
