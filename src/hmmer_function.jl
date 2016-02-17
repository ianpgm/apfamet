function run_hmmsearch(input_filenames, db_filename=Pkg.dir("apfamet")*"/db/Pfam-A.hmm", cores="1", sampleIDs=[], project_table_file = "project_table.txt", overwrite=false)
	#Check to see whether SampleIDs are present, if not then just number them starting from 1
	if sampleIDs == []
		sampleIDs = map(string,range(1,length(input_filenames)))
	end
	
	#If there are SampleIDs, check that there's enough of them
	if length(input_filenames) != length(sampleIDs)
		error("Number of Sample IDs and filenames do not match")
	end
	
	#Check to see whether hmmpress has been run, if not then run it
	if isfile(db_filename*".h3i") == false
		run(`hmmpress $db_filename`)
	end
	
	#Check to see whether project table file already exists, if it does and overwrite is set to the default false then stop the function
	if isfile(project_table_file)
		if overwrite==false
			error(project_table_file*" already exists, pick another name, delete the existing file, or specify overwrite=true")
		end
	end
	
	#Go through the input files and run hmmsearch, taking care not to overwrite any existing output files
	for filename in input_filenames
		output_filename = filename*".hmmsearch_tblout"
		if isfile(output_filename)
			if overwrite==false
				error(filename*".hmmsearch_tblout already exists, pick another name, delete the existing hmmsearch output, or specify overwrite=true")
			end
		end		
		run(`hmmsearch --tblout $output_filename --cpu $cores $db_filename $filename`)
	end
	
	#Write project table ready for the next step
	output_table = open(project_table_file, "w")
	write(output_table, "SampleID\tFilename\n")
	for (i, filename) in enumerate(input_filenames)
		write(output_table, sampleIDs[i]*"\t"*filename*".hmmsearch_tblout"*"\n")
	end
	
	close(output_table)
		
end
