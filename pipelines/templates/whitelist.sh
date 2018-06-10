umi_tools whitelist --stdin ${reads[0]} \\
--bc-pattern=${params.whitelist.bc_pattern} \\
--set-cell-number=${params.whitelist.cell_number} \\
--log2stderr > whitelist.txt;