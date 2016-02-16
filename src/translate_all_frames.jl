using Bio

function translate_all_frames(dna_fasta_files, format=Bio.Seq.FASTQ)
	for filename in dna_fasta_files
		input_file = open(filename, format)
		output_file = open(filename*".faa")
		for record in input_file
			print(record.name)
			write(output_file, ">"*record.name*".0"*"\n"*convert(AbstractString,Bio.Seq.translate(convert(Bio.Seq.RNASequence, record.seq), Bio.Seq.bacterial_plastid_genetic_code, false))*"\n")
		end
	end
end
		#output_file = open(filename*".faa")
		