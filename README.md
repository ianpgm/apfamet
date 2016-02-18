# apfamet
Apfamet is a [Julia](http://www.julialang.org/) package for analysing metagenomes based on the abundance of [HMMer](http://www.hmmer.org/) models (for example like those in [Pfam](http://pfam.xfam.org/)) in metagenomic datasets. It is intended for interactive use, for example from the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/) or [IJulia](https://github.com/JuliaLang/IJulia.jl) notebook. Apfamet should work on any platform Julia works on (Windows, Mac OS X, Linux), but has currently only been tested on Mac OS X. 

##Installation
1. Apfamet uses [HMMer](http://www.hmmer.org/) to search your sequences from HMM protein models - this software, specifically `hmmpress` and `hmmsearch`, needs to be installed and [in your PATH](http://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them) to work.
2. [Download and install the appropriate Julia version](http://julialang.org/downloads/) for your platform.
3. Open the [Julia REPL](http://docs.julialang.org/en/release-0.4/manual/interacting-with-julia/).
4. Type `Pkg.clone("https://github.com/ianpgm/apfamet", "apfamet")`
5. Before you can start using apfamet you will need a HMMer database. The default option is to use Pfam in the `apfamet/db` database. To find out how to download the most recent Pfam database and where exactly where to put this database file you can use the following commands in the REPL:

    ```
    using apfamet 
    apfamet.get_pfam_database()
    ```
6. To translate your sequences in all six-reading frames you will currently need to have [Python 2.X](https://www.python.org/) installed along with [Biopython](http://biopython.org/wiki/Main_Page). These dependencies will go away in future versions of apfamet.

##Tutorial
This flowchart provides an overview of apfamet's current functionality. Detailed instructions for using apfamet follow below.
![apfamet flowchart](https://github.com/ianpgm/apfamet/blob/master/doc/apfamet_flowchart.png)

1. Open the Julia REPL (or Jupyter notebook if you prefer)
2. Load the apfamet package you installed previously:
	
	```
	using apfamet
	```
3. Set your working directory to the location of your input (trimmed and dereplicated) sequences in fastq or fasta format and translate them in all six reading frames using the following command:
	
	```
	cd("/directory/with/my/sequencing_data/")
	apfamet.translate_sixframes(["metagenome_sample_1.fastq","metagenome_sample_2.fastq","metagenome_sample_3.fastq"])
	```
The default sequencing format is fastq (error encoding is irrelevant since it's ignored), if your sequences are in fasta format then you can specify that as follows:
	
	```
	apfamet.translate_sixframes(["metagenome_sample_1.fasta","metagenome_sample_2.fasta","metagenome_sample_3.fasta"], format="fasta")
	```
4. You will now have a protein fasta file for each of your input files with an extra `.faa` extension in the filename. Next comes the time-consuming step - HMMer. I imagine that many people might find it more convenient to run this step on a powerful server or cluster, but my instructions here will assume that you're doing everything on the same machine. To run your sequences against the default Pfam database, use the following command:
	
	```
	apfamet.run_hmmsearch(["metagenome_sample_1.fastq.faa","metagenome_sample_2.fastq.faa","metagenome_sample_3.fastq.faa"])
	```
There are a lot more optional variables that you can specify for `apfamet.run_hmmsearch`, this is the complete function:
	
	```
	apfamet.run_hmmsearch(input_filenames; db_filename=joinpath(Pkg.dir("apfamet"),"db","Pfam-A.hmm"), cores="1", sampleIDs=[], project_table_file = "apfamet_project_table.txt", overwrite=false)
	```
	
	+ `db_filename` the path of another hmm file if you don't want to use the default (most recent Pfam-A). Your replacement HMM database will have to have the necessary rpoB models from Pfam (see the section on normalisation below) for the calculation of rpoB equivalent values.
	+ `cores` to use more cores and make hmmsearch go faster
	+ `sampleIDs` to specify IDs for each of your sample (in the same order as the faa files for each of those samples). The default will be to number them starting from 1. You can also choose to change these sample IDs later
	+ `project_table_file` this is the input file for later steps - you will edit this to include metadata at a later point. The default filename is "apfamet_project_table.txt", but it will probably make more sense to choose a name that describes the collection of files you are processing.
	+ `overwrite` self-explanatory - if there's already output files with the names you've specified, should apfamet overwrite them? The default is not to overwrite them and generate an error.
5. You can now add metadata (for example, geochemical data) to your project table file - each parameter is a new column, with a name (no spaces) you give it in the header. This can be achieved by opening your project table file ("apfamet_project_table.txt" if you used the default) in a text editor (by putting tabs between the columns on each line) or Excel. [Click here](https://github.com/ianpgm/apfamet/blob/master/test/sample_metadata.txt) to see an example of this file. 
6. Now you need to read the hmmsearch results into Julia and carry out normalisation of the read counts. You need to find a temporary name for your project (`my_project` in this example) and run the following command (assuming default options were used in `run_hmmsearch()`):
	
	```
	my_project = apfamet.new_project()
	```
If you used a non-default HMM database or name for your project table file, these can be specified as follows:
	
	```
	my_project = apfamet.new_project(project_table_file="my_apfamet_project_table.txt", hmm_database="my_special_database.hmm")
	```
7. The first thing you should do once you make a new project is to save it to your hard drive. The following command will save four files to your working directory called "apfamet_project.apfamet_project_table", "apfamet_project.apfamet_read_counts_table", "apfamet_project.apfamet_hmm_database_info","apfamet_project.apfamet_rpoB_equiv_table". Each of these files are tab-delimited text files suitable for import into other programs (R, Excel, etc.).
	
	```
	apfamet.save_project(my_project)
	```
You can also specify a base filename other than the default and allow overwriting existing files with the same name as follows:
	
	```
	save_project(base_filename="my_apfamet_project", overwrite=true)
	```
To load a project that was previously saved use this command, specifying the `base_filename` prefix specified in the `save_project()` function:
	
	```
	my_project = load_project("my_apfamet_project")
	```
8. Now you have your project run through hmmsearch, loaded into memory and ready to analyse. To start out, why don't we try plotting the abundances of some Pfam models. The command below will plot the normalised (RpoB equivalent) abundance of the Pfam families of interest in each sample. The numbers at the end of each bar refer to the raw number of reads used to calculate the rpoB equivalent value. In this case, we produce the plot in PDF format - other available formats include PNG, SVG, and PS - see the [Gadfly documentation](http://gadflyjl.org/) for more guidance. The [Pfam website](http://pfam.xfam.org/) is the best resource for finding the IDs of relevant Pfam families.
	
	```
	myplot = apfamet.plotmodel(["MCR_alpha","Na_H_antiporter"], my_project)
	Gadfly.draw(Gadfly.PDF("myplot.pdf",5Gadfly.inch,5Gadfly.inch),myplot)
	```
![plotmodel()](https://github.com/ianpgm/apfamet/blob/master/doc/plotmodel_example.png)
9. How about plotting your HMM abundances against one of your metadata parameters? This is how you do that:
	
	```
	correlplot = apfamet.plotcorrelation(["MCR_alpha","Na_H_antiporter"],:CH4_mM,my_project)
	Gadfly.draw(Gadfly.PDF("correlplot.pdf",5Gadfly.inch,5Gadfly.inch),correlplot)
	```
![plotcorrelation()](https://github.com/ianpgm/apfamet/blob/master/doc/plotcorrelation_example.png)
10. Apfamet currently has a quick and crude principal-components-analysis method built in. It shows the two components with the most variation, but not how much variation those components actually show (like a scree plot). It's basically just good for getting a quick overview of your data. More correct multivariate statistical analysis can be carried out using Julia's [MultivariateStats package](https://github.com/JuliaStats/MultivariateStats.jl).
	
	```
	pca_plot = apfamet.perform_pca(my_project)
	Gadfly.draw(Gadfly.PDF("pca_plot.pdf",5Gadfly.inch,5Gadfly.inch),pca_plot)
	```
![perform_pca()](https://github.com/ianpgm/apfamet/blob/master/doc/perform_pca_example.png)

##How apfamet's normalisation works 
Apfamet uses the gene encoding RNA polymerase beta subunit (_rpoB_) to normalise counts in metagenomes. _rpoB_ is a universal, single-copy gene in prokaryotes, and thus a good basis for normalisation and comparison between data sets. There are seven pfam models apfamet uses to detect the RpoB protein sequence, these are:
```
["RNA_pol_Rpb2_1","RNA_pol_Rpb2_2","RNA_pol_Rpb2_3","RNA_pol_Rpb2_4","RNA_pol_Rpb2_45","RNA_pol_Rpb2_5","RNA_pol_Rpb2_6","RNA_pol_Rpb2_7"]
```
For each metagenomic sample, the number of reads for each of those models is found and "plotted" against the amino-acid length of each model, like so:

![normalisation example](https://github.com/ianpgm/apfamet/blob/master/doc/normalisation_plot.png)

Linear regression is performed (as represented by the line) using [Julia's GLM package](https://github.com/JuliaStats/GLM.jl) and the slope is taken as the value of reads/length for rpoB in that sample. Then, for every HMM found in the metagenomic sample, a reads/length figure is calculated based on the length of the HMM in the database and the number of reads that HMM hits in the metagenome. Then reads/length in the HMM is divided by reads/length in the rpoB, to define an abundance for each HMM in "rpoB equivalents". Since there's on average one rpoB in every genome, then one rpoB equivalent should be roughly equivalent to every genome in the metagenome having one copy of that HMM (on average). This unit (rpoB) makes metagenomes comparable based on a natural internal standard, in spite of different numbers of reads, qualities, and other variations from sample to sample.

##Roadmap
Features planned for apfamet in the future include:

+ Pearson correlation between model abundance and metadata variables.
+ A function to split projects based on metadata, sets of HMMs, and sets of sample IDs. Also a function to merge projects.
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
No, not at all! The tutorial in this file will tell you everything you need to know. It will help if you have some familiarity with interactive command-line based software like [R](https://www.r-project.org/), [mothur](http://www.mothur.org/) or [Qiime](http://qiime.org/), but my hope is that I've given adequate instructions for beginners - if not then please let me know. That being said some knowledge of Julia, especially its statistics packages and plotting capabilities, will enable you to do a lot more with apfamet's output than you would without that knowledge.

**Is your question not answered here? Don't hesitate to [submit an issue](https://github.com/ianpgm/apfamet/issues) or [contact me via email or telephone](http://pure.au.dk/portal/en/persons/id%2825504f0f-4132-4611-8697-0019cedc5d5d%29.html)**