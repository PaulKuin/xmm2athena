# XMM2athena-OM

**Introduction**

The codes found here were developed whilst enhancing the serendipitous UV source survey catalogue (SUSS) Version 5 of the XMM-Newton Optical Monitor by classifying the sources and investigating variability.  

The SUSS 5.0 catalogue is a catalogue where for each observed sky observation done at a given time the photometry in up to six filters is presented as well as related parameters like quality. This work made an intermediate product being the catalogue with just one record for each astronomical source (the single source catalogue). 

As nearby sources show parallax and proper motion, over time their detected postion will depend on the time (epoch) of observation. For this reason a subset of sources with proper motion larger than 10 milli-arcseconds per year was extracted from the Gaia catalogue. This means the positions of these sources can be found at any time using the proper motion and parallax in the Gaia catalogue. Each OM observation was matched to a subset of the Gaia catalogue for the epoch of that observation. This allows then to determine which records in the SUSS belong to a certain astronomical source for building the single source catalogue.

For each source a search was done for auxiliary data from ground-based surveys (PanStars,SkyMapper, UKIDS, VISTA and ALLWISE) which were used to augment the single source catalogue. Since for some sources multiple values of the magnitude in a certain band exist, these were used to find typical values (mean, error) to add for each source. 

For the classification subsets were determined of sources being a star, galaxy, or QSO. 

The classification method CLAXBOI from Hugo Tranin (https://github.com/htranin/classificationXray) was adapted for the UV-Optical data because CLAXBOI was written to classify the X-ray sources. Parameters used have been colours (magnitude differences) which are independent of distance. One magnitude difference between the Gaia-G magnitude and the ground-based g also depends on the extent of the source, since Gaia has a very small footprint on the sky, much smaller than groundbased observations, and turns out to be discriminating for galaxies. 

The choice of the method allowed to classify each source, trace which parameters affected the classification, and worked well with the very inhomogeneous catalogue: some sources have only a few magnitudes. Of course, the more magnitudes are present the better the classification. Also, sources fainter than the detection limit from Gaia-G (around > 20 magnitudes) were hampered to some degree in distinghuising QSO and galaxy. 

For variability, sources with multiple observations in a OM-SUSS band were selected. Only data with the best quality, quality flag = 0, were selected. Several other constraints have been applied to select good data.

**Code description**

The file named xmm2athena.py contains not only the python scripts, but also notes written while processing in steps. However, part of the processing was done using TopCat/stilts scripts. For those you will need to download Mark Taylors TopCat: www.star.bris.ac.uk/~mbt/topcat which is Java code that allows for very fast table operations and is preferred above, e.g., astropy.tables for speed on large catalogs. The scripts using stilts can be found because 'stilts' is in the name.  We found there is a bit of a learning curve for some of the operations, and having some familiarity with java can be useful. 

There are also some shell scripts. We used some of them for the matching of the observed sources in an observation to the Gaia DR3 sources using the on-line available service "Arches". See http://serendib2023.astro.unistra.fr/ARCHESWebService/XMatchARCHES.  While doing that for all the observations (OBSIDs) in SUSS 5.0, we had to restart several times so had to keep track of which OBSIDs had been completed. 

Python scripts to create the single source catalogue in two steps are present: single_source_cat.py with in the cats/ directory _v1 and _v2. These codes processed the SUSS sources with multiple observations by taking the median magnitude and associated error whilst selecting the quality 0 (best) data and determining the associated new quality flags.  The difference between v1 and v2 is mainly that the spread of input magnitudes was used to get a measure of varibility. In V2 that is done by computing the chi-squared, but using all data even of worse quality. The selection of only quality 0 data is done much later in xmm2athena variability processing. 
The second step of the processing is done after merging in the proper motion and parallax information by matching to the full Gaia DR3, including sources with hardly any proper motion, and also matching sources with a different ID but the same standard J2000, Epoch 2000 position. The code was also adapted so it can be used on the Swift UVOT catalogue (OM and UVOT were built the same way, with UVOT having different processing and Swift being in just a 90-minute orbit). The V1 files are obsolete but included for being how the earlier versions were processed. 

For more efficient processing of a large number of observations using arches is not optimal. The approach was changed to selecting all observations in a certain area, and matching using the stilts program. Here is also made use of the addEpoch2Source.sh script. 

In addition there are some shell scripts which were written for simple operations and some of those were not used in the processing at all.

**acknowledgement**
This project has received funding from the European Union's Horizon 2020 research and innovation programme under Grant Agreement No. 101004168. 



