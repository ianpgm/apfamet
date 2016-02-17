# apfamet
Apfamet is a [Julia](http://www.julialang.org/) package for analysing metagenomes based on the abundance of [HMMer](http://www.hmmer.org/) models (for example like those in [Pfam](http://pfam.xfam.org/)) in metagenomic datasets. It is intended for interactive use, for example from the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/) or [IJulia](https://github.com/JuliaLang/IJulia.jl) notebook. Apfamet should work on any platform Julia works on (Windows, Mac OS X, Linux), but has currently only been tested on Mac OS X. 

##Installation
1. [Download and install the appropriate Julia version](http://julialang.org/downloads/) for your platform.
2. Open the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/).
3. Type `Pkg.clone("https://github.com/ianpgm/apfamet", "apfamet")`
4. Before you can start using apfamet you will need a HMMer database and place its ".hmm" file in the `apfamet/db` directory. To find out how to download the most recent Pfam database and where exactly where to put this database file you can use the following commands in the REPL:
```
using apfamet
apfamet.get_pfam_database()
```
5. To translate your sequences in all six-reading frames you will currently need to have [Python 2.X](https://www.python.org/) installed along with [Biopython](http://biopython.org/wiki/Main_Page). These dependencies will go away in future versions of apfamet.

##Tutorial
This flowchart provides a basic overview of apfamet's current functionality. Detailed instructions for using apfamet follow below.
![apfamet flowchart](https://github.com/ianpgm/apfamet/blob/master/doc/apfamet_flowchart.png)

1. Open the Julia REPL (or Jupyter notebook if you prefer)
2. Load the apfamet package you installed previously:
```
using apfamet
```
3. Take your input (trimmed and dereplicated) sequences in fastq or fasta format and translate them in all six reading frames using the following command:
```
apfamet.translate_sixframes(["/path/to/file/my_metagenome_reads.fastq"])
```
You can also specify multiple files simultaneously like this:
```
apfamet.translate_sixframes(["/path/to/file/metagenome1.fastq","/path/to/file/metagenome2.fastq","/path/to/file/metagenome3.fastq"])
```
Or first set your working directory and then you don't need to specify the full path
```
cd("/path/to/file/")
apfamet.translate_sixframes(["metagenome1.fastq","metagenome2.fastq","metagenome3.fastq"])
```
4. You will now have a protein fasta file for each of your input files with an extra `.faa` extension in the filename. Next comes the time-consuming step - HMMer. For this step you will need to work directly from the command line (you can type `exit()` to exit Julia at any time) to run hmmsearch. I imagine that many people might find it more convenient to run this step on a shared server or cluster, but my instructions here will assume that you're doing everything on the same machine. You will also need to decide which HMM database to use at this point - apfamet was originally intended for use with Pfam, but should work with any HMM database. If you choose to use the default Pfam database that Julia comes with, you can find out where it is on your hard drive with the following command within the Julia REPL:
```
locate_default_database()
```
You first need to prepare that database using the `hmmpress` command:
```
hmmpress /path/to/database/database_file.hmm
```
...then search that database using your translated metagenomic reads as query sequences:
```
hmmsearch --tblout metagenome1.txt /path/to/database/database_file.hmm my_metagenome_reads.fastq.faa



##Roadmap
Features planned for apfamet in the future include:

+ Pearson correlation between model abundance and metadata variables.
+ Pairwise hypothesis testing to check for significant differences between groups of samples.
+ A search function for the HMM database, to search descriptions for keywords of interest
+ A function to export all reads with hits to a given model or set of models
+ Assembled data - checking for HMM abundance in assembled data is faster than checking every read, and data sets are getting too big for the read-based approach. That's why I would like to be able to use apfamet to combine assembled data with coverage information for analysis of assembled data.

Please [submit an issue](https://github.com/ianpgm/apfamet/issues) if you have suggestions for more features!

##Questions and Answers
**Why analyse metagenomes with HMMer instead of BLAST?**
The short answer is that there are already excellent tools out there for BLAST-based analysis, such as [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/). The long answer is that HMMer has some clear advantages over BLAST for analysing metagenomes. For example, for defining functions based on BLAST you need to set a sequence similarity (or BLAST score) cutoff. The trouble is that a sensible cutoff will depend on the function - some functions change with a single amino-acid substitution, sometimes two proteins with 30% sequence identity or less will share an identical function. HMMer distills a group of functionally related proteins down to its core essence (the hidden markov model - HMM) and bases its assessment on that, so that this varying degree of conservation between functions is taken into consideration. The other reason is that HMMer is generally better than BLAST for identifying distant homologues, which is usually the situation in metagenomes from environments with few cultured representatives.

**Why is apfamet written in [Julia](http://www.julialang.org/)?**
The functionality in this package was originally spread across some scripts written in R and Python, and when I discovered Julia it seemed like a convenient way to [combine the best of both worlds](http://schroed-mic.net/index.php/2015/12/17/thoughts-after-a-month-using-the-julia-programming-language/) and write the whole thing in one, easy-to-install, cross-platform package. Julia has a large and growing set of [statistical](http://juliastats.github.io/) and graphical tools that apfamet's output could be used with directly, so that seemed like a good way of giving users a lot of power in how they carried out their analysis.

**Do I need to be familiar with Julia to use apfamet?**
No, not at all! The tutorial in this file will tell you everything you need to know. It will help if you have some familiarity with command-line based tools like [mothur](http://www.mothur.org/) or [Qiime](http://qiime.org/), but my hope is that I've given adequate instructions for beginners - if not then please let me know. That being said some knowledge of Julia, especially its statistics packages and plotting capabilities, will enable you to do a lot more with apfamet's output than being without that knowledge.

**Is your question not answered here? Don't hesitate to [submit an issue](https://github.com/ianpgm/apfamet/issues) or [contact me via email or telephone](http://pure.au.dk/portal/en/persons/id%2825504f0f-4132-4611-8697-0019cedc5d5d%29.html)**