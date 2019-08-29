

#ifndef __SAC_H__
#define __SAC_H__

/** 
 * IIR Filter Prototypes 
 *
 * @see xapiir
 */
#define SAC_BUTTERWORTH       "BU"
#define SAC_BESSEL            "BE"
#define SAC_CHEBYSHEV_TYPE_I  "C1"
#define SAC_CHEBYSHEV_TYPE_II "C2"

/***
 * IIR Filter Types 
 * 
 * @see xapiir
 */
#define SAC_BANDPASS          "BP"
#define SAC_HIGHPASS          "HP"
#define SAC_LOWPASS           "LP"
#define SAC_BANDREJECT        "BR"

/** 
 * FIR Filter Type 
 *   
 * @see firtrn
 */
#define SAC_HILBERT           "HILBERT"
#define SAC_DERIVATIVE        "DERIVATIVE"

/** 
 * Window Type
 * 
 * @see crscor, window
 */ 
#define SAC_HAMMING           "HAMMING"
#define SAC_HANNING           "HANNING"
#define SAC_RECTANGLE         "RECTANGLE"
#define SAC_COSINE            "COSINE"
#define SAC_TRIANGULAR        "TRIANGULAR"

/** 
 * Filter
 *   
 *   Yes, this is what you want.  It is an IIR (Infinte Impulse Response) 
 *   filter and is the same that is used in lowpass, highpass, bandpass 
 *   and bandreject.  
 *
 *  @note Sensible Default
 *
 *  @bug This function needs to be wrapped in a filter function
 *       which is easier to use than what is presented here and
 *       functions in a similar manner to what is found in SAC. 
 *       Removal of the variables \p trbndw and \p a which are 
 *       only used with Cheyshev filters would be good, along with
 *       setting functions which only do lowpass, highpass, ....
 *        
 */
void xapiir ( float      data[], 
              int        nsamps, 
              char      *aproto, 
              double     trbndw, 
              double     a, 
              int        iord, 
              char      *type, 
              double     flo, 
              double     fhi, 
              double     ts, 
              int        passes);

/** 
 * Compute the envelope of a function
 *
 */
void envelope(int        n, 
              float     *in, 
              float     *out);

/** 
 * Finite Impulse Response Digital Filter 
 *
 * Can perform a hilbert transform or a derivative
 * 
 * @bug This function needs to be wrapped in a hilbert transform
 *      function to be made more accesible.  The same should be
 *      done for the derivative function as well.
 *
 */
void firtrn(char     *ftype, 
            float     x[], 
            int       n, 
            float     buffer[], 
            float     y[]);


/** 
 * Cross Correlation 
 * 
 *  Compute the cross correlation of two signals.
 *
 *  @note Sensible Defaults
 *     - \p nwin = 1
 *     - \p wlen = \p nsamps / \p nwin
 *     - \p type = SAC_RECTANGLE
 *     - \p nfft = Power of 2 greater than longest series
 *     
 *  @bug This function is capable of preforming a convolution
 *       if the input of either input data is reversed.
 *       This function needs to be wrapped in a more useable
 *       form for "convolution" and "correlation" with sensible
 *       defaults, as are set and used in the default SAC program
 */
void crscor(float     data1[], 
            float     data2[], 
            int       nsamps, 
            int       nwin, 
            int       wlen, 
            char     *type, 
            float     c[], 
            int      *nfft, 
            char     *err, 
            int       err_s);


/** 
 * Power of Two 
 *
 * Find the next largest power of two greater than \p num
 *
 */
int next2(int num);


#endif /* __SAC_H__ */
