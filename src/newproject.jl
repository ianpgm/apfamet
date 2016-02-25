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

function read_hmm_tabular_output(filename, ID,seqinfo)
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
			translation_frame = split(split_line[1], '.')[end]
			if haskey(matched_sequences_dict, target_sequence_name)
				if best_domain_evalue < matched_sequences_dict[target_sequence_name][1]
					matched_sequences_dict[target_sequence_name] = [best_domain_evalue,family,translation_frame]
				end
			else
				matched_sequences_dict[target_sequence_name] = [best_domain_evalue,family,translation_frame]
			end
		end
	end
	Pfam_dict = Dict()
	seq_df = DataFrame()
	for sequence in matched_sequences_dict
		if seqinfo
			newline_seq_df = DataFrame()
		end
		family = sequence[2][2]
		if haskey(Pfam_dict, family)
			Pfam_dict[family] = Pfam_dict[family] + 1
			if seqinfo
				newline_seq_df = DataFrame(seqID = sequence[1], translation = sequence[2][3], model = family, evalue = sequence[2][1])
			end
		else
			Pfam_dict[family] = 1
			if seqinfo
				newline_seq_df = DataFrame(seqID = sequence[1], translation = sequence[2][3], model = family, evalue = sequence[2][1])
			end
		end
		if seqinfo
			seq_df = vcat(seq_df, newline_seq_df)
		end
	end
	
	output_dict = Dict()
	output_dict["pfam_df"] = dict_to_df(Pfam_dict, ID)
	output_dict["seq_df"] = seq_df
	
	return output_dict
end

function add_newsample_to_pfam_table(filenames, IDs, pfam_table,seqinfo)
	new_pfam_table = pfam_table
	full_seq_df = DataFrame()
	for (i, filename) in enumerate(filenames)
		print("Reading "*filename*"\n")
		column_to_add = read_hmm_tabular_output(filename, IDs[i],seqinfo)
		new_pfam_table = join(column_to_add["pfam_df"], new_pfam_table, on = :PFAM_Model, kind = :outer)
		if seqinfo
			seq_df = column_to_add["seq_df"]
			seq_df[:Sample] = fill(IDs[i], size(seq_df,1))
			full_seq_df = vcat(full_seq_df, seq_df)
		end
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
	output_dict = Dict()
	output_dict["new_pfam_table"] = new_pfam_table
	output_dict["seq_df"] = full_seq_df
	return output_dict
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

#function import_hmm_output(project_table)
#	blank_table = DataFrame(PFAM_Model=[])
#	pfam_table = add_newsample_to_pfam_table(project_table[:Filename],project_table[:SampleID],blank_table)
#	output_dict = Dict()
#	output_dict["pfam_table"] = pfam_table["new_pfam_table"]
#	output_dict["seq_df"] = pfam_table[seq_df]
#	return output_dict
#end


function new_project(;project_table_file="apfamet_project_table.txt", hmm_database="", seqinfo=true)
	project_table = readtable(project_table_file, separator='\t')
	
	blank_table = DataFrame(PFAM_Model=[])
	pfam_table = add_newsample_to_pfam_table(project_table[:Filename],project_table[:SampleID],blank_table,seqinfo)
	read_counts_table = pfam_table["new_pfam_table"]
	seq_info_table = pfam_table["seq_df"]
	
	hmm_database_info = prepare_pfam_database(hmm_database)
	rpoB_equiv_table = reads_to_rpob_equiv(read_counts_table, hmm_database_info)
	return Dict("project_table"=>project_table,"read_counts_table"=>read_counts_table,"hmm_database_info"=>hmm_database_info,"rpoB_equiv_table"=>rpoB_equiv_table, "seq_info_table"=>seq_info_table)
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
	writetable(base_filename*".apfamet_seq_info_table", project["seq_info_table"], separator='\t')
end

function load_project(base_filename)
	project_table = readtable(base_filename*".apfamet_project_table", separator='\t')
	read_counts_table = readtable(base_filename*".apfamet_read_counts_table", separator='\t')
	hmm_database_info = readtable(base_filename*".apfamet_hmm_database_info", separator='\t')
	rpoB_equiv_table = readtable(base_filename*".apfamet_rpoB_equiv_table", separator='\t')
	seq_info_table = readtable(base_filename*".apfamet_seq_info_table", separator='\t')
	return Dict("project_table"=>project_table,"read_counts_table"=>read_counts_table,"hmm_database_info"=>hmm_database_info,"rpoB_equiv_table"=>rpoB_equiv_table, "seq_info_table"=>seq_info_table)
end