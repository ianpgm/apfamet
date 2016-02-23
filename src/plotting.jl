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
