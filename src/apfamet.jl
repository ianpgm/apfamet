module apfamet

using DataFrames
using GLM
using Gadfly
using Compose
using MultivariateStats
#using FTPClient
#using GZip

#Using a PyCall/Biopython solution to translation for the time being due to apparent bugs in BioJulia's translate function (see https://github.com/BioJulia/Bio.jl/issues/133)

include(Pkg.dir("apfamet")*"/src/pycall_translation_function.jl")
include(Pkg.dir("apfamet")*"/src/hmmer_function.jl")


function locate_default_database()
	print(joinpath(Pkg.dir("apfamet"),"db"))
end

function na_to_zero(input_value)
	if isna(input_value)
		return 0
	else
		return input_value
	end
end

function dict_to_df(input_dict, ID)
	new_df = DataFrame(PFAM_Model=collect(keys(input_dict)))
	new_df[symbol(ID)] = convert(Array{Int64},collect(values(input_dict)))
	return new_df
end

function read_hmm_tabular_output(filename, ID)
	hmm_output = open(filename)
	matched_sequences_dict = Dict()
	for line in eachline(hmm_output)
		if startswith(line, ['#'])
			continue
		else
			split_line = split(line)
			family = split_line[3]
			best_domain_evalue = float(split_line[8])
			target_sequence_name = join(split(split_line[1], '.')[1:end-1],'.')
			if haskey(matched_sequences_dict, target_sequence_name)
				if best_domain_evalue < matched_sequences_dict[target_sequence_name][1]
					matched_sequences_dict[target_sequence_name] = [best_domain_evalue,family]
				end
			else
				matched_sequences_dict[target_sequence_name] = [best_domain_evalue,family]
			end
		end
	end
	Pfam_dict = Dict()
	for sequence in matched_sequences_dict
		family = sequence[2][2]
		if haskey(Pfam_dict, family)
			Pfam_dict[family] = Pfam_dict[family] + 1
		else
			Pfam_dict[family] = 1
		end
	end
	
	return dict_to_df(Pfam_dict, ID)
end

function add_newsample_to_pfam_table(filenames, IDs, pfam_table)
	new_pfam_table = pfam_table
	for (i, filename) in enumerate(filenames)
		column_to_add = read_hmm_tabular_output(filename, IDs[i])
		new_pfam_table = join(column_to_add, new_pfam_table, on = :PFAM_Model, kind = :outer)
	end
#	for i in 1:size(new_pfam_table,1)
#		for j in 1:size(new_pfam_table,2)
#			if isna(new_pfam_table[i,j])
#				new_pfam_table[i,j] = 0
#			end
#		end
#	end
	for column in names(new_pfam_table)
		for (i, value) in enumerate(new_pfam_table[column])
			new_pfam_table[column][i] = na_to_zero(value)
		end
	end
	return new_pfam_table
end

function import_hmm_output(project_table)
	blank_table = DataFrame(PFAM_Model=[])
	pfam_table = add_newsample_to_pfam_table(project_table[:Filename],project_table[:SampleID],blank_table)
	return pfam_table
end

function dfsearch(df, colname, searchterm)
	hitlines = Int64[]
	for (i, value) in enumerate(df[colname])
		if contains(value,searchterm)
			push!(hitlines, i)
		end
	end
	return(df[hitlines,:])
end

function tabdict(filename)
	desc_table = Dict()
	for line in eachline(open(filename))
		desc_table[split(line, '\t')[1]] = strip(split(line, '\t')[2],'\n')
	end
	
	return desc_table
end


function reads_to_rpob_equiv(pfam_df, pfam_db_stats)
	#Get the names of the PFAM models for RpoB
	RpoB_models = ["RNA_pol_Rpb2_1","RNA_pol_Rpb2_2","RNA_pol_Rpb2_3","RNA_pol_Rpb2_4","RNA_pol_Rpb2_45","RNA_pol_Rpb2_5","RNA_pol_Rpb2_6","RNA_pol_Rpb2_7"]
	
	#Get their lengths
	RpoB_model_lengths = Array(Float64,0)
	for model in RpoB_models
		append!(RpoB_model_lengths, [float(pfam_db_stats[pfam_db_stats[:NAME] .== model, :LENG][1])])
	end
		
	#Calculate RpoB reads per aa for each sample
	RpoB_reads_per_aa = Array(Float64,0)
	for sample in names(pfam_df)
		RpoB_reads = Array(Float64,0)
		if sample != :PFAM_Model
			for RpoB_model in RpoB_models
				append!(RpoB_reads,float(pfam_df[pfam_df[:PFAM_Model] .== RpoB_model, sample]))
			end
			if length(RpoB_reads) != length(RpoB_model_lengths)
				print("It looks like there's some RpoB models missing from your dataset.")
			end
			
			
			data_for_linear_regression = DataFrame(reads=RpoB_reads,aa=RpoB_model_lengths)
			append!(RpoB_reads_per_aa,[coef(lm(reads~aa,data_for_linear_regression))[2]])
		end
	end
	
	#Calculate RpoB equivalents for every PFAM model in your dataset
	rpoB_table_output = DataFrame(PFAM_Model = [])
	for colnum in 2:size(pfam_df)[2]
		rpoB_table_output[names(pfam_df)[colnum]] = []
	end
	for rownum in 1:size(pfam_df)[1]
		model_name = pfam_df[rownum,:PFAM_Model]
		readnums = Array(pfam_df[rownum,2:size(pfam_df)[2]])
		if 	model_name in pfam_db_stats[:NAME]
			pfam_length = float(pfam_db_stats[pfam_db_stats[:NAME] .== model_name, :LENG][1])
		else
			error("At least one of the models in your hmmsearch output wasn't in the specified HMM database. Are you sure that the same HMM database was specified for both the run_hmmsearch() step and the new_project() step?\n")
		end
		reads_per_aa = []
		for readnum in readnums
			push!(reads_per_aa, readnum / pfam_length)
		end
		RpoB_equivalents = []
		for (i, reads) in enumerate(reads_per_aa)
			push!(RpoB_equivalents, reads / RpoB_reads_per_aa[i])
		end
		newrow = DataFrame(PFAM_Model = model_name)
		for colnum in 2:size(pfam_df)[2]
			newrow[names(pfam_df)[colnum]] = RpoB_equivalents[colnum-1]
		end
		rpoB_table_output = vcat(rpoB_table_output,newrow)
	end
	return rpoB_table_output
end

function new_project(;project_table_file="apfamet_project_table.txt", hmm_database="")
	project_table = readtable(project_table_file, separator='\t')
	read_counts_table = import_hmm_output(project_table)
	hmm_database_info = prepare_pfam_database(hmm_database)
	rpoB_equiv_table = reads_to_rpob_equiv(read_counts_table, hmm_database_info)
	return Dict("project_table"=>project_table,"read_counts_table"=>read_counts_table,"hmm_database_info"=>hmm_database_info,"rpoB_equiv_table"=>rpoB_equiv_table)
end

function save_project(project; base_filename="apfamet_project", overwrite=false)
	if isfile(base_filename)
		if overwrite == false
			error("Project "*base_filename*" already exists. Please select a new filename or set overwrite=true")
		else
			print("Overwriting previous project "*base_filename*" because overwrite=true.")
		end
	end
	writetable(base_filename*".apfamet_project_table", project["project_table"], separator='\t')
	writetable(base_filename*".apfamet_read_counts_table", project["read_counts_table"], separator='\t')
	writetable(base_filename*".apfamet_hmm_database_info", project["hmm_database_info"], separator='\t')
	writetable(base_filename*".apfamet_rpoB_equiv_table", project["rpoB_equiv_table"], separator='\t')
end

function load_project(base_filename)
	project_table = readtable(base_filename*".apfamet_project_table", separator='\t')
	read_counts_table = readtable(base_filename*".apfamet_read_counts_table", separator='\t')
	hmm_database_info = readtable(base_filename*".apfamet_hmm_database_info", separator='\t')
	rpoB_equiv_table = readtable(base_filename*".apfamet_rpoB_equiv_table", separator='\t')
	return Dict("project_table"=>project_table,"read_counts_table"=>read_counts_table,"hmm_database_info"=>hmm_database_info,"rpoB_equiv_table"=>rpoB_equiv_table)
end

function plotmodel(modelnames, project)
    data_to_plot = []
    key_labels =[]
	
	rpoB_eq_table = project["rpoB_equiv_table"]
	reads_table = project["read_counts_table"]
	hmm_database_info = project["hmm_database_info"]
	sample_names = map(string,names(project["rpoB_equiv_table"])[2:end])
	
    #Getting data out of the big table and putting it into a small data frame using the melt command to pivot that dataframe
    for (i, modelname) in enumerate(modelnames)
        single_model = melt(rpoB_eq_table[rpoB_eq_table[:PFAM_Model] .== modelname, :], :PFAM_Model)
        single_model_reads = melt(reads_table[reads_table[:PFAM_Model] .== modelname, :], :PFAM_Model)
		
        single_reads = []
        for readnum in single_model_reads[:value]
            push!(single_reads,string(readnum))
        end
        single_model[:reads] = single_reads
        single_model[:sample] = sample_names
		single_model[:PFAM_full] = hmm_database_info[hmm_database_info[:NAME] .== modelname, :DESC][1]
    relative_plot_position = 1:length(sample_names)
    relative_plot_position = relative_plot_position + 0.5 - i*(1/length(modelnames))
    single_model[:read_plot_pos] = relative_plot_position
        push!(data_to_plot, single_model)
    	push!(key_labels, hmm_database_info[hmm_database_info[:NAME] .== modelname, :DESC][1])
	end
    
    all_data_to_plot = data_to_plot[1]
    if length(data_to_plot) > 1
        for i in 2:length(data_to_plot)
            all_data_to_plot = vcat(all_data_to_plot, data_to_plot[i])
        end
    end
    
    #plotting
    modelplot = plot(all_data_to_plot, x=:value,y=:sample, colour=:PFAM_full, Geom.bar(position=:dodge,orientation=:horizontal),
    Guide.ylabel(""),
    Guide.xlabel("rpoB equivalents"),
    	Theme(bar_highlight=color(colorant"black"),
    	key_position=:bottom,
    	default_color=color(colorant"black"),
    	panel_stroke=color(colorant"black"),
    	grid_color=color(colorant"gray"),
    	major_label_font="Helvetica",
    	major_label_color=color(colorant"black"),
    	key_title_color=color(colorant"white"),
    	minor_label_font="Helvetica",
    	key_label_font="Helvetica",
    	minor_label_color=color(colorant"black"),
		point_label_font_size=6pt),
	Guide.annotation(compose(context(),
    text(all_data_to_plot[:value]+0.01,all_data_to_plot[:read_plot_pos],all_data_to_plot[:reads],[hleft])))
    )

	return modelplot
end

function perform_pca(project)
	
	rpoB_eq_table = project["rpoB_equiv_table"]
	sample_names = map(string,names(project["rpoB_equiv_table"])[2:end])
	
	rpoB_eq_matrix = convert(Array, rpoB_eq_table)
	rpoB_eq_matrix_converted = convert(Array{Float64,2},rpoB_eq_matrix[:,2:end])
	PCA_model = fit(PCA, rpoB_eq_matrix_converted)

	PC1 = transform(PCA_model, rpoB_eq_matrix_converted)[1,:]
	PC2 = transform(PCA_model, rpoB_eq_matrix_converted)[2,:]

	pca_plot = plot(x=PC1,y=PC2, Guide.annotation(compose(context(),text(PC1,PC2,map(string,sample_names),[hleft]))))
	
	return pca_plot
end

function plotcorrelation(modelnames,parameter,project)
    
	rpoB_eq_table = project["rpoB_equiv_table"]
	metadata_table = project["project_table"]
	hmm_database_info = project["hmm_database_info"]

	data_to_plot = []
    key_labels =[]

    #Getting data out of the big table and putting it into a small data frame using the melt command to pivot that dataframe
    for modelname in modelnames
        single_model = melt(rpoB_eq_table[rpoB_eq_table[:PFAM_Model] .== modelname, :], :PFAM_Model)
        single_model[:param] = reverse(metadata_table[parameter])
        single_model[:PFAM_full] = hmm_database_info[hmm_database_info[:NAME] .== modelname, :DESC][1]
        push!(data_to_plot, single_model)
	end

    #Putting all the data frames for each PFAM model together
    all_data_to_plot = data_to_plot[1]
    if length(data_to_plot) > 1
        for i in 2:length(data_to_plot)
            all_data_to_plot = vcat(all_data_to_plot, data_to_plot[i])
        end
    end
    
    #Plotting all of that data together
    plot(all_data_to_plot, x=:param,y=:value, color=:PFAM_full,
    Theme(default_color=color(colorant"black"),
    panel_stroke=color(colorant"black"),
    grid_color=color(colorant"gray"),
    major_label_font="Helvetica",
    major_label_color=color(colorant"black"),
    key_title_color=color(colorant"white"),
    minor_label_font="Helvetica",
    key_label_font="Helvetica",
    minor_label_color=color(colorant"black"),
    key_position=:bottom),
    Guide.xlabel(string(parameter)), 
    Guide.ylabel("rpoB equivalents")
    #Guide.annotation(compose(context(), text(rpoB_eq,metadata_subset,sample_names)))
    )
end



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
	
	print("This function is not yet working. Please download the file "*current_pfam_url*current_pfam_file*", unzip it and place it in "*target_location*".\n\nLinux and MacOS X users can use the following commands from the shell (not the Julia REPL - hit ';' to access the shell from within the Julia REPL):\ncurl "*current_pfam_url*current_pfam_file*" -o "*target_location*current_pfam_file*"\ngunzip "*target_location*current_pfam_file*"\n")
	
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

#Loading base data


end