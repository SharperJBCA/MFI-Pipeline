MFI-Pipeline
============

Collection of codes useful for processing QUIJOTE MFI data

Summary of modules:

The following are 'Toolbox' modules or Level-1 data-processing modules:
<ul>
<li> Binning: Contains a Fortran90 and Python wrapper routines for binning or down sampling data.</li>
<li> CalFitting: Contains routines for determining the average MFI calibration diode signal over an observations.</li>
<li> Coordinates: Fortran90 implementation of Denis' T-Point model and a Python wrapper to call it.</li>
<li> DataAccess: Some general Python routines to help with preparing file lists and also a routine for reading in MFI data files.</li>
<li> MFI_INFO: Contains key parameter files for calibrating MFI data. </li>
<li> Masking: Routines for masking emphemeris sources and noisey/spikey data. </li>
<li> SourceFitting: Routine for fitting Gaussian sources in TOD </li>
<li> WaveFitter: Routines for fitting periodic data that is discontinuous and not an integer number of wavelengths </li> 
</ul>

The Scipts module contains a number of Level-2 data-processing routines that combine one or more 'Toolbox' modules to complete a specific MFI data processing task.
