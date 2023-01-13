README This directory includes the observations used for the BTSI model as of November 8, 2022. Unless stated otherwise, the most recent version with most thorough post-processing offered by the instrument team is used. To create the observation matrix used to run the MCMC model, run makeobsmatrix.m from the main directory. Modifications, such as save path and demeaning, can be altered at the header of makeobsmatrix.m. Data files have been processed and homogenized to follow the following format:

Description of columns:
1st col: Julian day of observation
2nd col: Observed value in native unit of observer
format: (f7.1 f)




Original observations are downloaded from the following online sources:

NIMBUS/HF:

 "The Nimbus-7 Solar Total Irradiance:  A new Algorithm for its 
Deviation" by D.V. Hoyt, H.L. Kyle, J.R. Hickey, and R.H. Maschhoff (J.G.R.,
vol 97, pp.51-63)

https://www.ngdc.noaa.gov/stp/solar/solarirrad.html
Accessed: 2018/8/16

ACRIM1/SMM:

https://www.ngdc.noaa.gov/stp/solar/solarirrad.html
Accessed: 2018/8/16

ACRIM2/UARS:

https://www.ngdc.noaa.gov/stp/solar/solarirrad.html
Accessed: 2018/8/16

ACRIM3:

https://asdc.larc.nasa.gov/project/ACRIM%20III/ACR3L2DM_1
Accessed: 2018/8/17

ERBE/ERBS:

https://www.ngdc.noaa.gov/stp/solar/solarirrad.html
Accessed: 2018/8/17

PREMOS/PICARD:

http://idoc-picard.ias.u-psud.fr/sitools/client-user/res/html/projectdescription.html
Accessed: 2022/2/4

SORCE:
v19L3
http://lasp.colorado.edu/home/sorce/data/
Accessed: 2022/10/27

TCTE:
http://lasp.colorado.edu/home/tcte/data/
Accessed: 2021/12/8

TSIS-1:
http://lasp.colorado.edu/home/tsis/data/
Accessed: 2022/10/27

VIRGO/SOHO
The data are available from the ftp server ftp.pmodwrc.ch in the directory 
/data/irradiance/virgo/TSI/virgo_tsi_d_ and virgo_tsi_h_6_005_1805.dat
Accessed: 2022/10/27

SILSO:
https://wwwbis.sidc.be/silso/datafiles
Accessed: 2022/10/27

Bremen Mg-II:
Composite MgII index V5  mark.weber@uni-bremen.de
Accessed: 2022/10/27
