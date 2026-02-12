
rule combine_reports:
    input:
        reports = expand(rules.reformat_sumstat.output.report, locuseq=my_lb.locuseq)
    output:
        combined = ws_path("combined_report.tsv")
    resources:
        runtime=lambda wc, attempt: 60 + attempt * 30
    run:
        # Read each file and append to a list, then concatenate
        df_list = []
        for f in input.reports:
            # Check if file is not empty to avoid errors
            if os.path.getsize(f) > 0:
                temp_df = pd.read_csv(f, sep='\t')
                df_list.append(temp_df)
        
        # Merge all into one dataframe
        final_df = pd.concat(df_list, ignore_index=True)
        
        # Save to output
        final_df.to_csv(output.combined, sep='\t', index=False)
