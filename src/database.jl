function prepare_pfam_database(path_to_pfam_database="")
	if path_to_pfam_database==""
		path_to_pfam_database = joinpath(Pkg.dir("apfamet"),"db","Pfam-A.hmm")
	end
	pfam_stats_df = DataFrame(NAME=[],DESC=[],LENG=[])
	pfam_file = open(path_to_pfam_database)
	for line in eachline(pfam_file)
		if startswith(line, "NAME")
			global new_line = DataFrame(NAME=[split(line)[2]])
		elseif startswith(line, "DESC")
			new_line[:DESC] = strip(join(split(line)[2:end], ' '), '\n')
		elseif startswith(line, "LENG")
			new_line[:LENG] = parse(Int,split(line)[2])
			pfam_stats_df = vcat(pfam_stats_df, new_line)
		end
	end
	return pfam_stats_df
end

function get_pfam_database()
	current_pfam_url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
	current_pfam_file = "Pfam-A.hmm.gz"
	target_location = joinpath(Pkg.dir("apfamet"),"db")
	
	print("This function is not yet working. Please download the file "*current_pfam_url*current_pfam_file*", unzip it and place it in "*target_location*".\n\nLinux and MacOS X users can use the following commands from the shell (not the Julia REPL - hit ';' to access the shell from within the Julia REPL):\ncurl "*current_pfam_url*current_pfam_file*" -o "joinpath(target_location,current_pfam_file)"\ngunzip "joinpath(target_location,current_pfam_file)"\n")
	
#	print("Downloading current release of the Pfam database from "*current_pfam_url*"Pfam-A.hmm.gz\n")
	
#	ftp_init()
#	ftp_options = RequestOptions(url="ftp.ebi.ac.uk/pub/databases/Pfam/current_release/")
#	ftp_get("Pfam-A.hmm.gz", ftp_options, "Pfam-A.hmm.gz")
#	output_file = open(default_pfam_database_location*"Pfam-A.hmm-2.gz", "w")
#	write(output_file, read(ftp_stream.body))
	
#	ftp_object = ftp_get(host = "ftp.ebi.ac.uk/pub/databases/Pfam/current_release/")
#	binary(ftp_object)
#	download(ftp_object, "Pfam-A.hmm.gz",default_pfam_database_location*"Pfam-A.hmm.gz")
#	ftp_cleanup()
#	close(ftp_object)
	
	
#	print("Unzipping current release of the Pfam database to "*default_pfam_database_location*"Pfam-A.hmm.gz\n")
#	unzipped_output = open(default_pfam_database_location*"Pfam-A.hmm", "w")
#	gzipped_file = GZip.open(default_pfam_database_location*"Pfam-A.hmm.gz")
#	write(unzipped_output, readall(gzipped_file))
end
