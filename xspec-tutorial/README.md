# Using XSPEC for COSI spectral analysis and simulations

This tutorial is meant to serve as a starting point for using the X-ray spectral fitting package with COSI data. The main XSPEC website is https://heasarc.gsfc.nasa.gov/xanadu/xspec/ and a more extensive manual can be found there. XSPEC is installed by default as part of the HEASOFT tools. Thus, it is installed as part of the COSItools.  

XSPEC is an extremely mature software package that a large portion of the high-energy community uses. It has a huge number of models and flexibility for adding custom models as needed. Once the basic commands are known, it is easy to use, and spectral analysis can be performed efficiently. It has a very convenient spectral simulation tool that will be helpful to the COSI science team and will broaden the COSI community. Although not discussed in this tutorial, there is a pyXSPEC package that allows for automated scripts to be written.

Here, we assume that the user has extracted a spectral file (source_and_background.pha), a background file (background.pha), a response matrix (response.rmf), and an ancillary response file (area.arf) - these are provided in the above download. **Command#1** below gives an example of using the FTOOL `fdump` to see what is inside the file. The useful information is in the binary table extension. For example, there are two fields in each row of data (CHANNEL and COUNTS), and there are 502 rows in the table.  

The other thing that is important to check is that you have specified the names of the background file, the response matrix, and the ancillary response file in the header under the keywords BACKFILE, RESPFILE, and ANCRFILE.  When you load the file into XSPEC, these files will be read in automatically.  The easiest way to put the file names in the header is to use the FTOOL `fv`, search for the keywords, and type or copy the filenames in.

Both the `fdump` and `fv` commands are apart of FTOOLS and can be used in the regular command line or from inside XSPEC.

#### Command#1
```
[MacBook-Pro-4:cosi/2022/tutorial] jat% fdump source_and_background.pha
```
The output will be as follows:
```
Name of optional output file[STDOUT] 
Names of columns[-] 
Lists of rows[-] 
SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -32 / number of bits per data pixel
NAXIS   =                    0 / number of data axes
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
TELESCOP= 'COSI    '
INSTRUME= 'COSI    '
END
 
XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                    8 / width of table in bytes
NAXIS2  =                  502 / number of rows in table
PCOUNT  =                    0 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                    2 / number of fields in each row
TTYPE1  = 'CHANNEL '           / label for field   1
TFORM1  = 'J       '           / data format of field: 4-byte INTEGER
TTYPE2  = 'COUNTS  '           / label for field   2
TFORM2  = 'J       '           / data format of field: 4-byte INTEGER
TUNIT2  = 'count   '           / physical unit of field
EXTNAME = 'SPECTRUM'           / name of this binary table extension
TELESCOP= 'COSI    '           / telescope/mission name
INSTRUME= 'COSI    '           / instrument/detector name
EXTNAME = 'SPECTRUM'           / name of extension
FILTER  = 'NONE    '           / filter type if any
EXPOSURE=               38.912 / integration time in seconds
BACKFILE= 'background.pha'     / background filename
BACKSCAL=                   1. / background scaling factor
CORRFILE= 'NONE    '           / associated correction filename
CORRSCAL=                   1. / correction file scaling factor
RESPFILE= 'response.rmf'       / associated rmf filename
ANCRFILE= 'area.arf'           / associated arf filename
AREASCAL=                   1. / area scaling factor
STAT_ERR=                    0 / no statistical error specified
SYS_ERR =                    0 / no systematic error specified
GROUPING=                    0 / no grouping of the data has been defined
QUALITY =                    0 / no data quality information specified
HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'SPECTRUM'           / PHA dataset
HDUVERS = '1.2.1   '           / version of format
POISSERR=                    T / Poissonian errors to be assumed
CHANTYPE= 'PI      '           / channel type (PHA or PI)
DETCHANS=                  502 / total number of detector channels
HISTORY File modified by user 'jat' with fv  on 2022-06-06T22:37:16
END
 
        CHANNEL     COUNTS
                    count
      1           0           0
      2           1           0
      3           2           0
      4           3           0
      5           4           0
…
```

## Part 1: Standard spectral fitting with XSPEC

Add the following lines to your login file (.cshrc or similar).
```
setenv HEADAS /path/to/your/installation/heasoft-6.29/x86_64-apple-darwin19.6.0
alias heasoft "source $HEADAS/headas-init.csh"
```
Then, type `heasoft` to initialize HEASOFT and `xspec` to start an XSPEC session.
```
[MacBook-Pro-4:cosi/2022/tutorial] jat% xspec
```
You will see the xspec version returned. In our example, we have XSPEC 12.12.0:
```
                XSPEC version: 12.12.0
        Build Date/Time: Fri Sep  3 10:58:24 2021
```

Then, load the spectral file with the `data` command.
```
XSPEC12>data source_and_background.pha
```
The following information about the spectral file will be printed to the command line:
```
1 spectrum  in use
 
Spectral Data File: source_and_background.pha  Spectrum 1
Net count rate (cts/s) for Spectrum:1  1.590e+01 +/- 9.214e-01 (49.8 % total)
 Assigned to Data Group 1 and Plot Group 1
  Noticed Channels:  1-502
  Telescope: COSI Instrument: COSI  Channel Type: PI
  Exposure Time: 38.91 sec
 Using fit statistic: chi
 Using Background File                background.pha
  Background Exposure Time: 573 sec
 Using Response (RMF) File            response.rmf for Source 1
 Using Auxiliary Response (ARF) File  area.arf
```

