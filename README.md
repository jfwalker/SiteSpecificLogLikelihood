# SiteSpecificLogLikelihood
This Github contains the code necessary to run the MSWE and MGWE analyses for analyzing edges

Both programs are written in perl and should be relatively trivial to run, however, they rely
on phyx and raxml which need to be compiled and can be found here:

RAxML: ```https://github.com/stamatak/standard-RAxML```

https://www.ncbi.nlm.nih.gov/pubmed/24451623

Phyx: ```https://github.com/FePhyFoFum/phyx```

https://academic.oup.com/bioinformatics/article/33/12/1886/2975328

Both programs run by typing in the program name followed by the name of the config file
The vertebrate dataset is in the main folder to test it out. Each folder also contains
the code and the files used from the paper/results of the analyses. The programs with
examples are in the main folder. The data used for the examples is the Chiari et al. 2012
nucleotide data.

https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-10-65

**First step is to adjust the config file, for both the MSWE and MGWE the files are almost identical
apart from a few available options to run with the MSWE.**


1. The first settings are pxrmt,pxbp and pxrr and if it's installed to your bin you should
be able to leave it as is, if not you need to specify the path all the way to the name of
the program.

2. Next option is raxml, same as the phyx programs

3. Name of the outfile, the outfile will contain the info from the analysis, including
parameters, likelihood calcs, AIC, conflicts identified, gene counts etc... This
outfile will be printed to the main directory all is being run in.

4. For the species option, this is the edge/relationship you want to question. The species that make
up the clade you are questioning should all be comma separated. So in the example the edge you are testing
is the monophyly of: alligator,caiman,phrynops,caretta,chelonoidis_nigra,emys_orbicularis

5. The topologies option is asking how many topologies you want to test. The different relationships should
be in your file and be the first relationships in that file. The program will analyze them in that order. If you
have 3 you want to question, then put the 3 as your first three in you gene tree input file and put 3 as the topologies
to question. If you want to identify conflict and see if there are dominant topologies, then run it with zero and it
will try to identify your most commonly conflicting relationships.

6. The supermatrix is just the name of you supermatrix file, which must be fasta, if it's not you can use pxs2fa from the
phyx package to convert it prior to running the program.

7. This is you model file, it recognizes the length of genes based on DNA, genename = start-stop

8. The Set is your tree set, so all the trees you want to test with the topologies in question at the top

9. Thread is how many threads should be used for RAxML

10. The test option will run the program on 2 genes instead of the whole dataset

11. The verbose option prints out a bit more info, such as which trees had to have which species removed because
of missing data.

12. The Folder will be the name of an output folder for different temporary files the program generates (e.g. phyx logfile)

13. The secret mode can be used if you have run the program already and do not want it to recalculate the whole supermatrix
likelihood.

14. Other modes in the MSWE should probably be turned off, this spits out pretty much all info it can for the dataset, in
the case of the carnivory the site information is about 8Gb and that's only 13 taxa so it's probably better to leave it off


To run MGWE: 
```
perl MGWE_Calc.pl ExampleMGWE.config
```
To run MSWE
```
perl MSWE_Calc.pl ExampleMSWE.config
```


**Some Notes**

1. Depending on the version of RAxML, where it's compiled, order of your sequences etc... you use the likelihood values may be a bit different, this can have a minor influenceon site counts etc... It should not influence overall results unless your likelihood Calcs are super similar.

2. If sites have identical likelihoods for the topology in question they are not counted since they do not support either tree,
if the total sites counted seems lower than the total from the matrix, this could be why

3. The config file is very sensitive to being filled out correctly and it's important to check the outfile to make sure no
errors have been detected.

4. Outfile and the output folder should have different names from anything in your working directory or else they will be
deleted and written over

5. This repo contains the data from Walker et al. 2017 (http://www.amjbot.org/content/104/6/858.short) and Chiari et al. 2012 (link to paper above), as a result it is a pretty big repo since the exact analyses from the paper are in here to make sure everything is open.
