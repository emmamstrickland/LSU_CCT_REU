# Leveraging Graph Convolutional Networks to Predict Antibiotic Susceptibility to Bacterial Strains.

Inputs:
--------
- AMR_download: (str) title of downloaded text file with bacteria data from PATRIC database.
- anti_url: (str) url link to antibiotic data from STRING database to be mapped with PATRIC data by taxon ID.
- bact_proteins_link: (str) url link to bacteria protein data downloaded from STRING.
- sequences_link: (str) url link to protein sequence data from STRING to be used in GCN.

Example Usage: (main function)
------------------------------

AMR_download = 'PATRIC_genomes_AMR.txt'  
anti_url = 'https://stringdb-downloads.org/download/species.v12.0.txt'  
bact_proteins_link = 'https://stringdb-downloads.org/download/protein.links.v12.0/'  
sequences_link = 'https://stringdb-downloads.org/download/protein.sequences.v12.0/'

"""  
Mapping PATRIC susceptibility data to STRING data  
"""

antis = STRING_data(anti_url, ['#taxon_id', 'STRING_name_compact'])  
bact = AMR_data(AMR_download, ['genome_name','taxon_id', 'antibiotic', 'resistant_phenotype'])  
combined_data = mapping(bact, antis)  
make_csv(combined_data, 'mapped_data.csv')  
remove_incomplete_data('/work/estrick/main_project/mapped_data.csv', 'resistant_phenotype', '/work/estrick/main_project/clean_mapped_data.csv')  

"""  
Generating Graphs  
"""

bacteria_data(combined_data, bact_proteins_link)  
sequences(combined_data, sequences_link)  
unzip_files('bacteria_graphs')  
unzip_files('protein sequences')  
move_size_files('bacteria_graphs')  
move_size_files('protein_sequences')  
all_graphs('/work/estrick/main_project/bacteria_proteins', '/work/estrick/main_project/bacteria_graphs')

"""  
Creating histograms and other data statistics  
"""

susceptibility('/work/estrick/main_project/clean_mapped_data.csv')  
unique_instances('/work/estrick/main_project/clean_mapped_data.csv', 'taxon_id')  
stats = statistics('/work/estrick/main_project/bacteria_graphs', '/work/estrick/main_project/stats.csv')  
histograms('/work/estrick/main_project/stats.csv', '/work/estrick/main_project/histograms')  
graph_ratio('/work/estrick/main_project/bacteria_proteins', '/work/estrick/main_project/stats.csv', '/work/estrick/main_project/histograms')

print("All done!")

