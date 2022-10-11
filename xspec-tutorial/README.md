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
To see what the spectrum looks like, use the following commands:
`XSPEC12>iplot` (with this command, the use enters “pgplot”)
`PLT>device /xs` (this makes plots go to the computer screen; note that d /xs also works)
`PLT>time off` (turns off the timestamp at the bottom of the window)
`PLT>lw 5` (makes line widths five times thicker)
`PLT>lw 5 on 1` (makes the data group 1 point five times thicker)
`PLT>font roman` (uses times roman fonts rather than the default)
`PLT>plot`

You should see a spectrum that looks exactly like the one shown in Figure 1.
![Figure 1](Figures/source_and_background_spectrum.png) 

Now, you can leave pgplot by typing `exit or `quit`.   And for future reference, to end an XSPEC session and return to a regular command line, the command `exit` is used.

The next step is to bin the spectrum and to remove the channels outside the spectral range.  Here, we will use a standard set of bins (provided in groups8ch.dat).  Other techniques for binning are described in an Appendix. We bin the spectra and mark the bad channels using the FTOOL grppha (once again, this can be done from within XSPEC or from a regular command line):

```
XSPEC12>grppha source_and_background.pha
Please enter output filename[] source_and_background_grp8ch.pha

-------------------------
  MANDATORY KEYWORDS/VALUES
  -------------------------
  --------------------------------------------------------------------
  --------------------------------------------------------------------
  EXTNAME   - SPECTRUM        Name of this BINTABLE
  TELESCOP  - COSI            Mission/Satellite name
  INSTRUME  - COSI            Instrument/Detector
  FILTER    - NONE            Instrument filter in use
  EXPOSURE  - 38.912          Integration time (in secs) of PHA data
  AREASCAL  - 1.0000          Area scaling factor
  BACKSCAL  - 1.0000          Background scaling factor
  BACKFILE  - background.pha  Associated background file
  CORRSCAL  - 1.0000          Correlation scaling factor
  CORRFILE  - NONE            Associated correlation file
  RESPFILE  - response.rmf    Associated redistribution matrix file
  ANCRFILE  - area.arf        Associated ancillary response file
  POISSERR  - TRUE            Whether Poissonian errors apply
  CHANTYPE  - PI              Whether channels have been corrected
  TLMIN1    - 0               First legal Detector channel
  DETCHANS  - 502             No. of legal detector channels
  NCHAN     - 502             No. of detector channels in dataset
  PHAVERSN  - 1.2.1           OGIP FITS version number
  STAT_ERR  - FALSE           Statistical Error
  SYS_ERR   - FALSE           Fractional Systematic Error
  QUALITY   - TRUE            Quality Flag
  GROUPING  - FALSE           Grouping Flag
  --------------------------------------------------------------------
  --------------------------------------------------------------------
```
`GRPPHA[] bad 0-180` (these channels will not change as long as you are using groups8ch.dat)
`GRPPHA[] bad 377-501` (these channels will not change as long as you are using groups8ch.dat)
`GRPPHA[] group groups8ch.dat` 
`GRPPHA[] exit`
(Note that you will see text in the square brackets.  The text comes from the grppha parameter file. If the text in the square brackets is the command that you want, then you can press return instead of typing in the command.)

Now, read in the file that was created after grouping
`XSPEC12> data source_and_background_grp8ch.pha`
Remove the bins outside the spectral range with
`XSPEC12> ignore bad`
We also convert the x-axis to energy using
`XSPEC12> setplot energy`

We plot the spectrum using
`XSPEC12> plot ldata`
`XSPEC12> iplot`
(and the other commands in pgplot above) and Figure 2 shows the result.

![Figure 2](Figures/source_and_background_grp8ch_spectrum.png)


