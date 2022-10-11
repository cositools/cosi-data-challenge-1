#Using XSPEC for COSI spectral analysis and simulations

This tutorial is meant to serve as a starting point for using the X-ray spectral fitting package with COSI data. The main XSPEC website is https://heasarc.gsfc.nasa.gov/xanadu/xspec/ and a more extensive manual can be found there. XSPEC is installed by default as part of the HEASOFT tools. Thus, it is installed as part of the COSItools.  

XSPEC is an extremely mature software package that a large portion of the high-energy community uses. It has a huge number of models and flexibility for adding custom models as needed. Once the basic commands are known, it is easy to use, and spectral analysis can be performed efficiently. It has a very convenient spectral simulation tool that will be helpful to the COSI science team and will broaden the COSI community. Although not discussed in this tutorial, there is a pyXSPEC package that allows for automated scripts to be written.

Here, we assume that the user has extracted a spectral file (source_and_background.pha), a background file (background.pha), a response matrix (response.rmf), and an ancillary response file (area.arf) - these are provided in the above download. Command#1 below gives an example of using the FTOOL “fdump” to see what is inside the file. The useful information is in the binary table extension. For example, there are two fields in each row of data (CHANNEL and COUNTS), and there are 502 rows in the table.  

The other thing that is important to check is that you have specified the names of the background file, the response matrix, and the ancillary response file in the header under the keywords BACKFILE, RESPFILE, and ANCRFILE.  When you load the file into XSPEC, these files will be read in automatically.  The easiest way to put the file names in the header is to use the FTOOL “fv,” search for the keywords, and type or copy the filenames in.

Both the fdump and fv commands are apart of FTOOLS and can be used in the regular command line or from inside XSPEC.

