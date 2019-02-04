# WRF4.0.1_Thompson
WRF V4.0.1 model runs using the Thompson-Eidhammer microphysics scheme for the Microphysics of Antarctic Clouds campaign. Simulations to test the influence of CCN number concentrations on cloud microphysics, using the case study observations from 27 November 2015 for comparison.

Issues with v3.9.1.1 -- error with mismatch between p_top in namelist and in input file. Not obvious why this error occurred.

------------------------------------------------

Success with v4.0.1 on ws8. (3-Jan-2019)

Trialing on Archer. (3-Jan-2019)

              -- Option 50 (dmpar: INTEL (ftn/icc): Cray XC) 
              
Successful compilation and runs without segmentation faults when chem/ code compiled (7-Jan-2019)
  
-- Model grinds to a halt at d02 2015-11-26_00:00:00 calling inc/HALO_EM_SBM_inline.inc. Does not abort with error message, continues to run. No change in filesize over ~30minutes, so definitely not doing anything (might have been an issue with the rsl.error output). Likely a memory error

4-Feb-2019: Decide to revert back to standard version of the model (without Nisg80 diagnostics). No graupel or snow number concentrations output from Thompson-Eidhammer scheme as default, so no validation checks possible without the extra work of outputting these. For ease (and to allow study to be conducted on ARCHER), revert back to default.
