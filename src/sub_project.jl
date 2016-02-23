function sub_project_by_sample(project, sample_names)
	columns = append!([:PFAM_Model], reverse(map(symbol, sample_names)))
	new_read_counts_table = project["read_counts_table"][columns]
	new_rpoB_equiv_table = project["rpoB_equiv_table"][columns]
	new_project_table = DataFrame()
	for sample in sample_names
		new_project_table = vcat(new_project_table, project["project_table"][project["project_table"][:SampleID] .== sample,:])
	end
	new_seq_info_table = DataFrame()
	for sample in sample_names
		new_seq_info_table = vcat(new_seq_info_table, project["seq_info_table"][project["seq_info_table"][:Sample] .== sample,:])
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
		new_seq_info_table = vcat(new_seq_info_table, project["seq_info_table"][project["seq_info_table"][:model] .== model,:])
	end
	
	return Dict("project_table"=>project["project_table"],"read_counts_table"=>new_read_counts_table,"hmm_database_info"=>project["hmm_database_info"],"rpoB_equiv_table"=>new_rpoB_equiv_table, "seq_info_table"=>new_seq_info_table)
	
	
end

function sub_project_by_metadata(project, parameter, expression)

end

function merge_projects(projects)
	
end	