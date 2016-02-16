# apfamet
Apfamet is a [Julia](http://www.julialang.org/) package for analysing metagenomes based on the abundance of [HMMer](http://www.hmmer.org/) models (for example like those in [Pfam](http://pfam.xfam.org/)) in metagenomic datasets. It is intended for interactive use, for example from the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/) or [IJulia](https://github.com/JuliaLang/IJulia.jl) notebook. Apfamet should work on any platform Julia works on (Windows, Mac OS X, Linux), but has currently only been tested on Mac OS X.

##Installation
1. [Download and install the appropriate Julia version](http://julialang.org/downloads/) for your platform.
2. Open the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/).
3. Type `Pkg.clone("https://github.com/ianpgm/apfamet", "apfamet")`
4. Before you can start using apfamet you will need a HMMer database and place its ".hmm" file in the `apfamet/db` directory. To find out exactly where to put this database file you can use the following commands in the REPL:
```
using apfamet
apfamet.get_pfam_database()
```

##Tutorial

