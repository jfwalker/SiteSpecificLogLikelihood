# SiteSpecificLogLikelihood
This Github contains the code necessary to run the MSWE and MGWE analyses for analyzing edges

Both programs are written in perl and should be relatively trivial to run, however, they rely
on phyx and raxml which need to be compiled and can be found here:

RAxML: ```https://github.com/stamatak/standard-RAxML```

Phyx: ```https://github.com/FePhyFoFum/phyx```


Both programs run by typing in the program name followed by the name of the config file
The vertebrate dataset is in the main folder to test it out. Each folder also contains
the code and the files used from the paper/results of the analyses.

First step is to adjust the config file

The first setting is pxrmt, and their just write the path to pxrmt


and run: 
```
perl MGWE_Calc.pl ExampleMGWE.config
```
