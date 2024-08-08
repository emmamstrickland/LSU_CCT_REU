# Leveraging Graph Convolutional Networks to Predict Antibiotic Susceptibility to Bacterial Strains.

Inputs:
--------
- AMR_download: (str) title of downloaded text file with bacteria data from PATRIC database.
- anti_url: (str) url link to antibiotic data from STRING database to be mapped with PATRIC data by taxon ID.
- bact_proteins_link: (str) url link to bacteria protein data downloaded from STRING.
- sequences_link: (str) url link to protein sequence data from STRING to be used in GCN.

Example Usage: (in main function)
------------------------------

    AMR_download = 'PATRIC_genomes_AMR.txt'  
    anti_url = 'https://stringdb-downloads.org/download/species.v12.0.txt'  
    bact_proteins_link = 'https://stringdb-downloads.org/download/protein.links.v12.0/'  
    sequences_link = 'https://stringdb-downloads.org/download/protein.sequences.v12.0/'  
    
    """  
    Gathering data from PATRIC and STRING  
    """
    
    antis = STRING_data(anti_url, ['#taxon_id', 'STRING_name_compact'])  
    bact = AMR_data(AMR_download, ['genome_name','taxon_id', 'antibiotic', 'resistant_phenotype'])  
    combined_data = mapping(bact, antis)  
    make_csv(combined_data, 'name_for_mapped_data_file.csv')  
    remove_incomplete_data('/path/to/unedited/mapped_data.csv', 'resistant_phenotype', '/path/to/edited/mapped_data.csv')

    """  
    Gathering files for GCN and creating graphs  
    """
    
    bacteria_data(combined_data, bact_proteins_link, 'folder_name_to_save_bacteria_downloads')  
    protein_sequences = sequences(combined_data, sequences_link, 'folder_name_to_save_protein_sequence_downloads')  
    unzip_files('folder_with_bacteria_downloads')  
    unzip_files('folder_with_protein_sequence_downloads')  
    move_size_files('folder_with_bacteria_size_files')  
    move_size_files('folder_with_sequences_size_files')  
    all_graphs('/folder/with/bacteria_protein_downloads', '/folder/to/save/bacteria_graphs')

    """  
    Creating histograms and other data statistics  
    """
    
    susceptibility('/path/to/edited_mapped_data.csv')  
    unique_instances('/path/to/edited_mapped_data.csv', 'taxon_id')  
    stats = statistics('path/to/folder/with/bacteria_graphs', '/path/to/save/graph_stats.csv')  
    histograms('/path/to/graph_stats.csv', '/path/to/save/histograms')  
    graph_ratio('/folder/with/bacteria_protein_downloads', '/path/to/graph_stats.csv', '/path/to/save/histograms')

    
    print("All done!")  
