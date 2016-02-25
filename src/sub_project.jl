function sub_project_by_sample(project, sample_names)
	columns = append!([:PFAM_Model], reverse(map(symbol, sample_names)))
	new_read_counts_table = project["read_counts_table"][columns]
	new_rpoB_equiv_table = project["rpoB_equiv_table"][columns]
	new_project_table = DataFrame()
	for sample in sample_names
		new_project_table = vcat(new_project_table, project["project_table"][project["project_table"][:SampleID] .== sample,:])
	end
	new_seq_info_table = DataFrame()
	if size(project["seq_info_table"]) != (0,0)
		for sample in sample_names
			new_seq_info_table = vcat(new_seq_info_table, project["seq_info_table"][project["seq_info_table"][:Sample] .== sample,:])
		end
	else
		new_seq_info_table = project["seq_info_table"]
	end
	
	return Dict("project_table"=>new_project_table,"read_counts_table"=>new_read_counts_table,"hmm_database_info"=>project["hmm_database_info"],"rpoB_equiv_table"=>new_rpoB_equiv_table, "seq_info_table"=>new_seq_info_table)
	
end

function sub_project_by_model(project, model_names)
	new_read_counts_table = DataFrame()
	new_rpoB_equiv_table = DataFrame()
	new_seq_info_table = DataFrame()
	for model in model_names
		new_read_counts_table = vcat(new_read_counts_table, project["read_counts_table"][project["read_counts_table"][:PFAM_Model] .== model,:])
		new_rpoB_equiv_table = vcat(new_rpoB_equiv_table, project["rpoB_equiv_table"][project["rpoB_equiv_table"][:PFAM_Model] .== model,:])
		if size(project["seq_info_table"]) != (0,0)
			new_seq_info_table = vcat(new_seq_info_table, project["seq_info_table"][project["seq_info_table"][:model] .== model,:])
		else
			new_seq_info_table = DataFrame()
		end
	end
	
	return Dict("project_table"=>project["project_table"],"read_counts_table"=>new_read_counts_table,"hmm_database_info"=>project["hmm_database_info"],"rpoB_equiv_table"=>new_rpoB_equiv_table, "seq_info_table"=>new_seq_info_table)
	
	
end

function sub_project_by_metadata(project, parameter, operator, value)
	expression = Expr(:call,operator,project["project_table"][parameter],value)
	new_sample_IDs = project["project_table"][eval(expression),:][:SampleID]
	return sub_project_by_sample(project, new_sample_IDs)
end

function merge_projects(project1, project2)
	sample_intersect = intersect(project1["project_table"][:SampleID],project2["project_table"][:SampleID])
	if length(sample_intersect) > 0
		error("Both projects share the following sample names and therefore cannot be merged:\n"*sample_intersect)
	end
	
	if project1["hmm_database_info"] != project2["hmm_database_info"]
		error("It looks like these projects were generated using different HMM databases. Merging them probably isn't the best idea.")
	end
	
	new_read_counts_table = join(project1["read_counts_table"], project2["read_counts_table"], on = :PFAM_Model, kind = :outer)
	new_rpoB_equiv_table = join(project1["rpoB_equiv_table"], project2["rpoB_equiv_table"], on = :PFAM_Model, kind = :outer)
	new_seq_info_table = vcat(project1["seq_info_table"],project2["seq_info_table"])
	new_project_table = vcat(project1["project_table"],project2["project_table"])
	new_hmm_database_info = project1["hmm_database_info"]
	
	return Dict("project_table"=>new_project_table,"read_counts_table"=>new_read_counts_table,"hmm_database_info"=>new_hmm_database_info,"rpoB_equiv_table"=>new_rpoB_equiv_table, "seq_info_table"=>new_seq_info_table)
	
end	