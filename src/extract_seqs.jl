using Bio

function extract_seqs(seq_IDs, filenames, format)
	output_array = []
	for filename in filenames
		print("Checking "*filename)
		for entry in open(filename, format)
			print(record.id)
			if record.id in seq_IDs
				push!(output_array, record)
			end
		end
	end
	return output_array
end

function extract_nt_seqs(models,project,filenames,output;format=Bio.Seq.FASTQ)
	seqs_to_retrieve = []
	for model in models
		append!(seqs_to_retrieve, project["seq_info_table"][project["seq_info_table"][:model] .== model, :seqID])
	end
	output_sequences = extract_seqs(seqs_to_retrieve, filenames, format)
	output =open(output, format)
	write(output, output_sequences)
end

function extract_aa_seqs(models,project,filenames,output;format=Bio.Seq.FASTA)
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
	output =open(output, format)
	write(output, output_sequences)
end