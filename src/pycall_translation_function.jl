using PyCall

@pyimport __builtin__
@pyimport Bio.SeqIO as SeqIO
@pyimport Bio.Seq as Seq
@pyimport Bio.Alphabet as IUPAC

function translate_sixframes(filenames, format="fastq")
	for filename in filenames
		file = open(filename)
		output_file = open(filename*".faa", "w")

		for record in SeqIO.parse(file, format)
			nt_sequence = record[:seq]
			sequence_id = record[:id]
			
			write(output_file, ">"*sequence_id*".0\n"*__builtin__.str(pycall(nt_sequence[:translate],PyAny, table=11))*"\n")
			
			frame_1_str = __builtin__.str(nt_sequence)[2:end]
			frame_1_seq = Seq.Seq(frame_1_str, IUPAC.generic_dna)
			write(output_file, ">"*sequence_id*".1\n"*__builtin__.str(pycall(frame_1_seq[:translate],PyAny, table=11))*"\n")
			
			frame_2_str = __builtin__.str(nt_sequence)[3:end]
			frame_2_seq = Seq.Seq(frame_2_str, IUPAC.generic_dna)
			write(output_file, ">"*sequence_id*".2\n"*__builtin__.str(pycall(frame_2_seq[:translate],PyAny, table=11))*"\n")
			
			rc_nt_sequence = pycall(nt_sequence[:reverse_complement], PyAny)
			write(output_file, ">"*sequence_id*".0rc\n"*__builtin__.str(pycall(rc_nt_sequence[:translate],PyAny, table=11))*"\n")
			
			frame_1rc_str = __builtin__.str(rc_nt_sequence)[2:end]
			frame_1rc_seq = Seq.Seq(frame_1rc_str, IUPAC.generic_dna)
			write(output_file, ">"*sequence_id*".1rc\n"*__builtin__.str(pycall(frame_1rc_seq[:translate],PyAny, table=11))*"\n")
			
			frame_2rc_str = __builtin__.str(rc_nt_sequence)[3:end]
			frame_2rc_seq = Seq.Seq(frame_2rc_str, IUPAC.generic_dna)
			write(output_file, ">"*sequence_id*".2rc\n"*__builtin__.str(pycall(frame_2rc_seq[:translate],PyAny, table=11))*"\n")			
			
		end

	end
end