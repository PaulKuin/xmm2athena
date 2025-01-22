# XMM2athena-OM
The codes found here were developed whilst enhancing the serendipitous UV source survey catalogue (SUSS) Version 5 of the XMM-Newton Optical Monitor by classifying the sources and investigating variability.  

The SUSS 5.0 catalogue is a catalogue where for each observed sky observation done at a given time the photometry in up to six filters is presented as well as related parameters like quality. This work made an intermediate product being the catalogue with just one record for each astronomical source (the single source catalogue). 

As nearby sources show parallax and proper motion, over time their detected postion will depend on the time (epoch) of observation. For this reason a subset of sources with proper motion larger than 10 milli-arcseconds per year was extracted from the Gaia catalogue. This means the positions of these sources can be found at any time using the proper motion and parallax in the Gaia catalogue. Each OM observation was matched to a subset of the Gaia catalogue for the epoch of that observation. This allows then to determine which records in the SUSS belong to a certain astronomical source for building the single source catalogue.

For each source a search was done for auxiliary data from ground-based surveys. 
