using PyCall

@pyimport __builtin__
@pyimport Bio.SeqIO as SeqIO
@pyimport Bio.Seq as Seq

function extract_seqs(seq_IDs, filenames, format)
	output_array = []
	for filename in filenames
		print("Checking "*filename*"\n")
		for record in SeqIO.parse(filename, format)
			if record[:id] in seq_IDs
				push!(output_array, record)
			end
		end
	end
	return output_array
end

function extract_nt_seqs(models,project,filenames,output;format="fastq")
	if size(project["seq_info_table"]) == (0,0)
		error("This project contains no sequence information table. Rerun new_project() with seqinfo=true to be able to extract sequences")
	end
	
	seqs_to_retrieve = []
	for model in models
		append!(seqs_to_retrieve, project["seq_info_table"][project["seq_info_table"][:model] .== model, :seqID])
	end
	print(seqs_to_retrieve)
	output_sequences = extract_seqs(seqs_to_retrieve, filenames, format)
	output_handle = __builtin__.open(output, "w")
	SeqIO.write(output_sequences,output_handle,format)
end

function extract_aa_seqs(models,project,filenames,output;format="fasta")
	if size(project["seq_info_table"]) == (0,0)
		error("This project contains no sequence information table. Rerun new_project() with seqinfo=true to be able to extract sequences")
	end
	
	nt_seqs_to_retrieve = []
	translations = []
	aa_seqs_to_retrieve = []
	for model in models
		append!(nt_seqs_to_retrieve, project["seq_info_table"][project["seq_info_table"][:model] .== model, :seqID])
		append!(translations, project["seq_info_table"][project["seq_info_table"][:model] .== model, :translation])
		for (i, seqID) in enumerate(nt_seqs_to_retrieve)
			push!(aa_seqs_to_retrieve, join([seqID,translations[i]],'.'))
		end
	end
	output_sequences = extract_seqs(aa_seqs_to_retrieve, filenames, format)
	output_handle = __builtin__.open(output, "w")
	SeqIO.write(output_sequences,output_handle,format)
end