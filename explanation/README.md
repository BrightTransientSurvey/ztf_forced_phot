# ZTF fps

The IPAC forced photometry service (fps) produces fixed-position, "forced" photometry on existing difference images in the ZTF IPAC archive. The [IPAC fps documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf) includes all the relevant details, here we provide a brief summary of the fps before describing the BTS light curve processing procedure.

 The fps takes a user-supplied position (and date range) and queries the ZTF subtraction image database for all images that contain the specified coordinates. PSF photometry is then measured at the specified position using the already-determined PSF model for the subtraction image. Results from these measurements are then returned in a standard IPAC table. There are many columns in the table, some of which are useful for flagging unreliable flux measurements, but the most important columns are: `jd`, `filter`, `zpdiff`, `forcediffimflux`, and `forcediffimfluxunc`.

While the fps provides both a zero-point and flux measurement, potential systematic offsets in the gain of the reference and new science images requires correction prior to measuring a calibrated flux. 

# BTS light curves

IPAC fps output requires a "baseline correction" ([see the docs](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf)). This baseline is best determined using a large (>~ 25) statistical sample of subtraction images that do not include emission from the transient of interest. Most transients have a suitable set of "pre-transient" images, but a significant fraction (e.g., many ZTF discoveries from 2018) require images taken long after the transient has faded. A suitable baseline can be identified visually, though with thousands of BTS transients we automate this procedure.

### Time of maximum

For individual transients we determine the time of maximum by calculating a running median with a 14 day window for every field + filter + chip + quadrant combination (designated `fcqfid` see [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5)) in the fps products. The final time of maximum is determined by taking the mean time of maximum for any g-band or r-band fcqfids taken on the primary ZTF observing grid (see [Bellm et al. (2018)](http://dx.doi.org/10.1088/1538-3873/aaecbe) for more details on the primary and secondary ZTF observing grids).

### Baseline correction

By default, the baseline for each `fcqfid` is determined using observations obtained more than *100 d prior* to the time of maximum light. This threshold is extremely conservative as very few transients have a rise time >100 d. If there are not enough suitable pre-transient observations (i.e., < 25) to measure the baseline correction, the baseline is determined using observations obtained more than *500 d after* the time of maximum light. Again this choice is conservative, though there are transients that have long rise times or that can be detected for years after they are discovered. The baseline correction `C` is calculated via the median `forcediffimflux` value in the baseline region.

We report the total number of observations used to calculate `C` and whether the baseline is measured in observations before or after maximum light. We advise caution when using light curves with few baseline observations (< ~25), and do not provide baseline corrections when there are < 10 observations in the baseline region. If there are > 25 pre-maximum baseline observations or (> 10 pre-maximum AND < 25 post-maximum baseline observations then we use the pre-maximum baseline. Otherwise we use the post-maximum baseline. Pre-maximum baseline measurements are always better than post-maximum measurements because there is always some transient emission present at late times (emission that is statistically insignificant after 500 d, hopefully). 

### Scaling the uncertainties

Following the recommendation in the [fps documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf), we re-scale the uncertainties of all observations based on the chi-squared distribution of observations in the baseline region. In brief, flux measurements in an "empty" patch of sky should follow a gaussian distribution and "excess" scatter points to some systematic uncertainty that has not been properly quantified (subtractions near galaxy nuclei often exhibit such excess scatter). If the chi-squared distribution is >1, then the uncertainties *of all observations* are multiplied by the square root of the reduced chi squared in the baseline region (see also [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5)).

### Flagging unreliable observations

We attempt to flag observations that are not reliable, in case users wish to remove these from their analysis. For the most part, such images were taken in extremely poor observing conditions and it was not possible to properly calibrate the PSF model. In particular, we flag images that match the following criteria as `poor_condutions` within the data products: 

1. `infobits` > 0
2. `scisigpix` > 25
3. `sciinpseeing` > 5

(see [Yao et al. (2019)](http://dx.doi.org/10.3847/1538-4357/ab4cf5) for further details). A non-zero value for `infobits` typically indicates a problem with the IPAC data processing pipeline. `infobits` > 33554432 corresponds to "cloudy" data (see Section 2.4 of the [ZTF Science Data System(ZSDS) Advisories & Cautionary Notes](http://web.ipac.caltech.edu/staff/fmasci/ztf/extended_cautionary_notes.pdf)). 

Note that observations that have been flagged as unreliable are excluded from any baseline calculations.

### Calibrated fluxes

To produce the final calibrated flux measurements, the baseline value `C` is subtracted from all flux measurements `forcediffimflux`. The uncertainties are scaled (see above), and the flux in microJy is calculated as:

10^(29 - 48.6/2.5 - 0.4 * `zpdiff` * (`forcediffimflux` - `C`))

## Final product

The final output from this procedure includes a csv file `ZTFYYnnnnnnn_fnu.csv` with nine columns:

1. JD - julian data of the observation
2. fnu_microJy - the flux in microJy
3. fnu_microJy_unc - the uncertainty on the flux in microJy
4. passband - the filter of the observations
5. fcqfid - field + filter + ccd + quadrant combination (see above)
6. passband - the filter (ZTF_g, ZTF_r, or ZTF_i)
7. N_baseline - the number of baseline observations used to calculate `C`
8. pre_or_post - =-1 if the baseline is calculated pre-maximum; =1 if the baseline is calculated post-maximum
9. poor_conditions - =1 if the observation is unreliable (see above)