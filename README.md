# ztf_forced_phot
Python routines to generate flux calibrated photometry from the ZTF forced photometry service provided by IPAC


## ZTF Forced Photometry

The IPAC data processing team for the [Zwicky Transient Facility (ZTF)](https://www.ztf.caltech.edu/) recently released a [forced photometry service (fps)](https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi) that provides fixed-position PSF photometry on all publicly available ZTF images. Details of the fps and the resulting outputs are summarized in the [documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf). 

This repository provides the software used to convert the raw output from the IPAC fps to calibrated flux measurements for all the transients found by the [Bright Transient Survey (BTS)](https://sites.astro.caltech.edu/ztf/bts/bts.php) during the first phase of ZTF (roughly 01 May 2018 –– 31 Oct 2020). The [explanation/](explanation/) folder includes some visual representations of the processing procedure, though we also refer the interested user to the IPAC [documentation](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf). The final light curves from this procedure are published in Miller et al. 2021, in prep.










#### Acknowledgements

The ZTF forced-photometry service was funded under the Heising-Simons Foundation grant #12540303 (PI: Graham).