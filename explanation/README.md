# ZTF fps

The IPAC forced photometry service (fps) produces fixed-position, "forced" photometry on existing difference images in the ZTF IPAC archive. The [IPAC fps documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf) includes all the relevant details, here we provide a brief summary of the fps before describing the BTS light curve processing procedure.

 The fps takes a user-supplied position (and date range) and queries the ZTF subtraction image database for all images that contain the specified coordinates. PSF photometry is then measured at the specified position using the already-determined PSF model for the subtraction image. Results from these measurements are then returned in a standard IPAC table. There are many columns in the table, some of which are useful for flagging unreliable flux measurements, but the most important columns are: `jd`, `filter`, `zpdiff`, `forcediffimflux`, and `forcediffimfluxunc`.

While the fps provides both a zero-point and flux measurement, potential systematic offsets in the gain of the reference and new science images requires correction prior to measuring a calibrated flux. 

# BTS light curves

IPAC fps output requires a "baseline correction" ([see the docs](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf)). This baseline is best determined using a large (>~ 25) statistical sample of subtraction images that do not include emission from the transient of interest. Most transients have a suitable set of "pre-transient" images, but a significant fraction (e.g., many ZTF discoveries from 2018) require images taken long after the transient has faded. A suitable baseline can be identified visually, though with thousands of BTS transients we automate this procedure.

### Flagging unreliable observations

We attempt to flag observations that are not reliable, in case users wish to remove these from their analysis. For the most part, such images were taken in extremely poor observing conditions and it was not possible to properly calibrate the PSF model. In particular, we flag images that match the following criteria as `poor_conditions` within the data products: 

1. `infobits` > 0
2. `scisigpix` > 25
3. `sciinpseeing` > 5

(see [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5) for further details). A non-zero value for `infobits` typically indicates a problem with the IPAC data processing pipeline. `infobits` > 33554432 corresponds to "cloudy" data (see Section 2.4 of the [ZTF Science Data System (ZSDS) Advisories & Cautionary Notes](http://web.ipac.caltech.edu/staff/fmasci/ztf/extended_cautionary_notes.pdf)). 

Note that observations that have been flagged as unreliable are excluded from any baseline calculations.

<img src="./../images/flagged_obs.jpg" raw=True>

### Time of maximum

For individual transients we determine the time of maximum by calculating a running median with a 14 day window for every field + filter + chip + quadrant combination (designated `fcqfid` see [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5)) in the fps products. The final time of maximum is determined by taking the mean time of maximum for any g-band or r-band fcqfids taken on the primary ZTF observing grid (see [Bellm et al. (2018)](http://dx.doi.org/10.1088/1538-3873/aaecbe) for more details on the primary and secondary ZTF observing grids).

### Baseline correction
By default, the baseline for each `fcqfid` is determined using observations obtained more than *100 d prior* to the time of maximum light and long after the transient has faded below the ZTF detection threshold. This time is (conservatively) determined by assuming the peak luminosity of the transient is purely powered by radioactive $^{56}$Co, and we then calculate the time to fade to 22.5 mag in the observed frame assuming the transient is at z = 0.09 (this choice of redshift is conservative). The baseline correction `C` is calculated via the median `forcediffimflux` value in the baseline region, and uncertainty of the baseline `C` is obtained using bootstrapping.

We report the total number of observations used to calculate `C` and whether the baseline is measured in observations before, after, or before and after the maximum light. We advise caution when using light curves with few baseline observations (< ~25), and do not provide baseline corrections when there are < 2 observations in the baseline region. If there are pre-maximum light or post-maximum light emission in the baseline, we provide "long rise" or "long decline" warning respectively. If both pre- and post- maximum light emissions are present, a "bad baseline" warning is given.

The figure below again shows SN 2019yvq, the running median used to determine the time of maximum and the baseline region used to determine `C`.

<img src="./../images/baseline_max.jpg" raw=True>

### Scaling the uncertainties

Following the recommendation in the [fps documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf), we re-scale the uncertainties of all observations based on the chi-squared distribution of observations in the baseline region. In brief, flux measurements in an "empty" patch of sky should follow a Gaussian distribution and "excess" scatter points to some systematic uncertainty that has not been properly quantified (subtractions near galaxy nuclei often exhibit such excess scatter). If the mean of the chi-squared distribution is < 1.5, then the uncertainties *of all observations* are multiplied by the square root of the median chi-squared distribution in the baseline region. Otherwise, if there are "good observations" (`infobits` = 0, see above), a [SuperSmoother](https://github.com/jakevdp/supersmoother/) model is used to obtain a smooth fit of chi-squared distribution as function of flux of these "good observations". The chi-squared distribution for the "bad observations" (`infobits` > 0) are obtained by interpolating the chi-squared distributions of the "good observations". In case the model fails, 1000s rolling medians on the "good observations" are used for the interpolation. The uncertainties *of all observations* are then multiplied by the square root of all the interpolated chi-squared distributions (see also [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5)).

<img src="./../images/scale_uncertainties.jpg" raw=True>

### Calibrated fluxes

To produce the final calibrated flux measurements, the baseline value `C` is
subtracted from all flux measurements `forcediffimflux`. The uncertainties are
scaled (see above), and the flux in microJy is calculated as:

<img src="https://render.githubusercontent.com/render/math?math=\Large f_\nu = 10^{29 - 48.6/2.5 - 0.4*\mathrm{zpdiff}}*(\mathrm{forcediffimflux} - C)">

## Final product

The final output from this procedure includes a csv file `ZTFYYnnnnnnn_fnu.csv` with nine columns:

1. JD - Julian date of the observation
2. fnu_microJy - the flux in microJy
3. fnu_microJy_unc - the uncertainty on the flux in microJy
4. passband - the filter of the observations (ZTF_g, ZTF_r, or ZTF_i)
5. programid - Program identifier (0 = engineering; 1 = public; 2 = private; 3 = Caltech time) ([see IPAC fps documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf))
6. fcqfid - field + filter + ccd + quadrant combination (see above)
7. zpdiff - photometric zeropoint for difference image
8. sys_unc_factor - systematic uncertainty used for scaling (see above)
9. poor_conditions - =1 if the observation is unreliable (see above)