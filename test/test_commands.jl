project_table = DataFrames.readtable("sample_metadata.txt", separator='\t')
pfam_table = apfamet.import_hmm_output(project_table)
pfam_stats_df = apfamet.prepare_pfam_database()
test_rpoB_equiv_table = apfamet.reads_to_rpob_equiv(pfam_table, pfam_stats_df)