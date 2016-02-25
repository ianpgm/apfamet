using HypothesisTests

function pearson_correlation(project,modelname,parameter)
     rpoB_eq_table = project["rpoB_equiv_table"]
	metadata_table = project["project_table"]
	single_model = melt(rpoB_eq_table[rpoB_eq_table[:PFAM_Model] .== modelname, :], :PFAM_Model)
     single_model[:param] = reverse(metadata_table[parameter])
    return cor(single_model[:param],single_model[:value])
end

function mann_whitney_compare_projects(project1,project2,modelname)
	project_1_vector = vec(Array{Float64}(project1["rpoB_equiv_table"][project1["rpoB_equiv_table"][:PFAM_Model] .== modelname,:][2:end]))
	project_2_vector = vec(Array{Float64}(project2["rpoB_equiv_table"][project2["rpoB_equiv_table"][:PFAM_Model] .== modelname,:][2:end]))
	return MannWhitneyUTest(project_1_vector, project_2_vector)
end

function mann_whitney_compare_models(project,model1,model2)
	model_1_vector = vec(Array{Float64}(project["rpoB_equiv_table"][project["rpoB_equiv_table"][:PFAM_Model] .== model1,:][2:end]))
	model_2_vector = vec(Array{Float64}(project["rpoB_equiv_table"][project["rpoB_equiv_table"][:PFAM_Model] .== model2,:][2:end]))
	return MannWhitneyUTest(model_1_vector, model_2_vector)
end
