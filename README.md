# MatchLocate2
Detect and locate small events from continuous seismic waveforms using templates (MatchLocate2)

                    MatchLocate2.0  2019/08/29
             https://github.com/Dal-mzhang/MatchLocate2
                Miao Zhang, Dalhousie Univesrity
                    E-mail: miao.zhang@dal.ca
                    
# Usage
(type "MatchLocate2", step by step tutorials can be found in examples)

Usage: MatchLocate2 -F(refevla/refevlo/refevdp) -R(maxlat/maxlon/maxh) -I(dlat/dlon/dh)
       -T(window/before/after) -H(CC/N(*MAD)) -D(INTD) -B(low/high) -O(ouput) INPUT.in
       
-F: searching center (e.g., 37.799/139.998/7.8).
-R: searching area (e.g., 0.05/0.05/5.0).
-I: searching interval (e.g., 0.01/0.01/1.0).
-T: time length of the reference phase (e.g., 4.0/1.0/3.0).
-H: cross-correlation thresholds CC && NMAD (e.g., 0.3/10.0, or 0.3/0.0, or 0.0/10.0).
-D: keep one event within INTD sec (e.g., 6.0).
-B: bandpass fitering for both templates and traces (e.g., 2/8).
-O: output (1,2,3) or don't output (0) the cross-correlogram or CC coefficient.
INPUT.in: directories of templates and continuous data, horizontal and vertical slowness, etc.


# Introduction:
Compared to the current methods of small event detection (template matching/matched filter), the Match&Locate places event detection to a lower magnitude level and extends the capability of detecting small events that have large distance separations from the template. The method has little dependence on the accuracy of the velocity models used, and, at the same time, provides high-precision location information of the detected small-magnitude events.

# References:
Zhang M. and Wen L. An effective method for small event detection: match and locate (M&L). Geophysical Journal International, 200 (3), 1523-1537, 2015, https://doi.org/10.1093/gji/ggu466

Zhang M. and Wen L. Seismological Evidence for a Low‐Yield Nuclear Test on 12 May 2010 in North Korea. Seismological Research Letters, 86 (1), 138-145, 2015, https://doi.org/10.1785/02201401170

Zhang M. and Wen L. Earthquake characteristics before eruptions of Japan's Ontake volcano in 2007 and 2014. Geophysical Research Letters, 42 (17), 6982–6988, 2015, https://doi.org/10.1002/2015GL065165

Our GPU version is available at https://github.com/MinLiu19/GPU-MatchLocate1.0; See reference: Liu M., Li H., Zhang M. and Wang T. Graphics Processing Unit-based Match and Locate (GPU-M&L): an improved Match and Locate method and its application. Seismological Research Letters, 91(2A), 1019-1029, 2020, https://doi.org/10.1785/0220190241 

A recent applicaiton of high-precision earthquake detection and location using M&L; See reference: Liu M., Zhang M. and Li H. Using the match-and-locate method to characterize the foreshocks of the July 2019 Mw 6.4 Ridgecrest, CA earthquake, Earthquake Science, 2022, https://doi.org/10.1016/j.eqs.2022.08.005
