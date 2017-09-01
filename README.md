# peterplots
Simple script for visualizing genome assembly quality
Usage in terminal from directory containing fasta files: 

	'Rscript --vanilla peter_plot_shovill.R *.fasta'

Will create output folder with renamed assemblies, as well as a pdf and interactive html containing plots showing length, coverage and gc content of each assembly's contigs. 

