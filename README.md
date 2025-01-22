# XMM2athena-OM
The codes found here were developed whilst enhancing the serendipitous UV source survey catalogue (SUSS) Version 5 of the XMM-Newton Optical Monitor by classifying the sources and investigating variability.  

The SUSS 5.0 catalogue is a catalogue where for each observed sky observation done at a given time the photometry in up to six filters is presented as well as related parameters like quality. This work made an intermediate product being the catalogue with just one record for each astronomical source (the single source catalogue). 

As nearby sources show parallax and proper motion, over time their detected postion will depend on the time (epoch) of observation. For this reason a subset of sources with proper motion larger than 10 milli-arcseconds per year was extracted from the Gaia catalogue. This means the positions of these sources can be found at any time using the proper motion and parallax in the Gaia catalogue. Each OM observation was matched to a subset of the Gaia catalogue for the epoch of that observation. This allows then to determine which records in the SUSS belong to a certain astronomical source for building the single source catalogue.

For each source a search was done for auxiliary data from ground-based surveys (PanStars,SkyMapper, UKIDS, VISTA and ALLWISE) which were used to augment the single source catalogue. Since for some sources multiple values of the magnitude in a certain band exist, these were used to find typical values (mean, error) to add for each source. 

For the classification subsets were determined of sources being a star, galaxy, or QSO. 

The classification method CLAXBOI from Hugo Tranin was adapted for the UV-Optical data because CLAXBOI was written to classify the X-ray sources. Parameters used have been colours (magnitude differences) which are independent of distance. One magnitude difference between the Gaia-G magnitude and the ground-based g also depends on the extent of the source, since Gaia has a very small footprint on the sky, much smaller than groundbased observations, and turns out to be discriminating for galaxies. 

The choice of the method allowed to classify each source, trace which parameters affected the classification, and worked well with the very inhomogeneous catalogue: some sources have only a few magnitudes. Of course, the more magnitudes are present the better the classification. Also, sources fainter than the detection limit from Gaia-G (around > 20 magnitudes) were hampered to some degree in distinghuising QSO and galaxy. 

For variability, sources with multiple observations in a OM-SUSS band were selected. Only data with the best quality, quality flag = 0, were selected. Several other constraints have been applied to select good data.



