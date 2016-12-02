# Readme for Terahertz Correlation Spectroscopy
Eric J. Rees, University of Cambridge. 2016. CC-BY-4.0

EJR Orcid ID: orcid.org/0000-0001-7478-5961

Contents
1. Published paper
2. Sample data
3. Matlab code for sample data analysis
4. Notes on future development


## 1. Published Paper
   Terahertz correlation spectroscopy infers particle velocity and rheological properties
   Eric J. Rees, Ke Su, and Axel Zeitler

   Opt. Lett. 41(14), 3289-3292 (2016)

   Optics Letters: https://www.osapublishing.org/ol/abstract.cfm?uri=ol-41-14-3289
   doi: 10.1364/OL.41.003289

   Openaccess.cam identifier and url: OA-8992
   https://www.repository.cam.ac.uk/handle/1810/256542

## 2. Sample data:
   Sample data for the THz Correlation Spectroscopy paper is archived at:

   data.cam: https://www.repository.cam.ac.uk/handle/1810/254844

   This data (~600 MB) is needed to run the script that generates the published figures

## 3. Matlab code for data analysis.

(a) THz_Corr_Vel_Code_Figs.m

This script generates all the figures that present data in the paper.

    This script was developed in Matlab 2013b.
    It reads in raw sample data from the opendata folder.
    To use this script, open Matlab, set the opendata folder as the working
    directory, and run this script.

(b) TOF.m

Time-of-flight analysis

    Analyses viscometer data to produce time-and-depth resolved
    velocimetry results using the terahertz correlation spectroscopy
    method presented in the paper.


(c) NOT INCLUDED IN THIS FOLDER - but useful for making greyscale compatible figures in the paper:

    cubehelix.m
    Implement's Dave Green's cubehelix colormap.
    Used for Figure 4.
    Cubehelix is a colormap with maximised constrast that can be
    converted to grayscale without loss of information.

    The cubehelix colormap is fully presented in:

    Green, D. A., 2011, `A colour scheme for the display of astronomical intensity images', Bulletin of the Astronomical Society of India, 39, 289.
    (2011BASI...39..289G at ADS. http://adsabs.harvard.edu/abs/2011BASI...39..289G )

    http://astron-soc.in/bulletin/11June/289392011.pdf

## 4. Notes on future development - tbd
