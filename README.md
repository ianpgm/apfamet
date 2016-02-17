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

##Tutorial
This flowchart provides a basic overview of apfamet's current functionality. Detailed instructions for using apfamet follow below.
[[https://github.com/ianpgm/apfamet/blob/master/doc/apfamet_flowchart.png]]


##Roadmap
Features planned for apfamet in the future include:
+Pearson correlation between model abundance and metadata.
+Pairwise hypothesis testing to check for significant differences between groups of samples.
+A search function for the HMM database, to search descriptions for keywords of interest
+A function to export all reads with hits to a given model or set of models
+Assembled data - checking for HMM abundance in assembled data is faster than checking every read, and data sets are getting too big for the read-based approach. That's why I would like to be able to use apfamet to combine assembled data with coverage information for analysis of assembled data.
Please [submit an issue](https://github.com/ianpgm/apfamet/issues) if you have suggestions for more features!

##Questions and Answers
**Why analyse metagenomes with HMMer instead of BLAST?**
The short answer is that there are already excellent tools out there for BLAST-based analysis, such as [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/). The long answer is that HMMer has some clear advantages over BLAST for analysing metagenomes. For example, for defining functions based on BLAST you need to set a sequence similarity (or BLAST score) cutoff. The trouble is that a sensible cutoff will depend on the function - some functions change with a single amino-acid substitution, sometimes two proteins with 30% sequence identity or less will share an identical function. HMMer distills a group of functionally related proteins down to its core essence (the hidden markov model) and bases its assessment on that, so that this varying degree of conservation between functions is taken into consideration. The other reason is that HMMer is generally better than BLAST for identifying distant homologues, which is usually the situation in metagenomes from environments with few cultured representatives.
**Why is apfamet written in [Julia](http://www.julialang.org/)?**
The functionality in this package was originally spread across some scripts written in R and Python, and when I discovered Julia it seemed like a convenient way to [combine the best of both worlds](http://schroed-mic.net/index.php/2015/12/17/thoughts-after-a-month-using-the-julia-programming-language/) and write the whole thing in one, easy-to-install, cross-platform package. Julia has a large and growing set of [statistical](http://juliastats.github.io/) and graphical tools that apfamet's output could be used with directly, so that seemed like a good way of giving users a lot of power in how they carried out their analysis.
**Is your question not answered here? Don't hesitate to [submit an issue](https://github.com/ianpgm/apfamet/issues) or [contact me via email or telephone](http://pure.au.dk/portal/en/persons/id%2825504f0f-4132-4611-8697-0019cedc5d5d%29.html)**