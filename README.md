# SiteSpecificLogLikelihood
This Github contains the code necessary to run the MSWE and MGWE analyses for analyzing edges

Both programs are written in perl and should be relatively trivial to run, however, they rely
on phyx and raxml which need to be compiled and can be found here:

RAxML: ```https://github.com/stamatak/standard-RAxML```

Phyx: ```https://github.com/FePhyFoFum/phyx```


Both programs run by typing in the program name followed by the name of the config file
The vertebrate dataset is in the main folder to test it out. Each folder also contains
the code and the files used from the paper/results of the analyses.

First step is to adjust the config file, for both the MSWE and MGWE the files are almost identical
apart from a few available options to run with the MSWE.


The first settings are pxrmt,pxbp and pxrr and if it's installed to your bin you should
be able to leave it as is, if not you need to specify the path all the way to the name of
the program.

Next option is raxml, same as the phyx programs

Name of the outfile, the outfile will contain the info from the analysis, including
parameters, likelihood calcs, AIC, conflicts identified, gene counts etc... This
outfile will be printed to the main directory all is being run in.




To run MGWE: 
```
perl MGWE_Calc.pl ExampleMGWE.config
```
To run MSWE
```
perl MSWE_Calc.pl ExampleMSWE.config
```
