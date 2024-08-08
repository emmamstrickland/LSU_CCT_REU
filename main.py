"""
Author:        Emma Strickland
PI:            Dr. Michal Brylinski
Title:         main.py
Date created:  06/13/24
Date updated:  08/07/24
"""
import requests
import csv
from bs4 import BeautifulSoup
import os
import gzip
import shutil
import numpy as np
import pandas as pd
from collections import defaultdict
from statistics import mean
import networkx as nx
import matplotlib.pyplot as plt
from networkx.utils.decorators import not_implemented_for
from concurrent.futures import ThreadPoolExecutor, as_completed

def AMR_data(file_name, included_columns):
    """
    Creates list of dictionaries from downloaded file from PATRIC database.
    
    file_name: (str) name of downloaded txt file in the same dictionary.
    included_columns: (str) list of the columns in txt file to be included in the dictionaries.
    """
    with open(file_name, 'r') as file:
        headers = file.readline().strip().split('\t')

        included_indices = [headers.index(col) for col in included_columns if col in headers]

        total_bact = 0
        dict_list = []
        print("Gathering bacteria data from PATRIC...")
        for line in file:
            total_bact += 1
            values = line.strip().split('\t')
            row_dict = {headers[i]: values[i] if i < len(values) else None for i in included_indices}
            dict_list.append(row_dict)
        print('Total instances from PATRIC:', total_bact)
#        print('PATRIC instances with needed columns', len(dict_list))
    return dict_list

def STRING_data(url, included_columns=None):
    """
    creates list of dictionaries from STRING downloads.
    
    included_columns: (str) optional list of columns in downloaded data to be included in list
    of dictionaries. If included_columns is None, function includes all columns
    """
    response = requests.get(url)

    if response.status_code == 200:
        content = response.text
        lines = content.splitlines()
        headers = lines[0].strip().split('\t')

        if included_columns is not None:
            included_columns = set(included_columns)

        dict_list = []
        print("Gathering antibiotic data from STRING...")
        for line in lines[1:]:
            values = line.strip().split('\t')
            row_dict = {headers[i]: values[i] for i in range(len(headers)) if included_columns is None or headers[i] in included_columns}
            dict_list.append(row_dict)
    else:
        raise Exception(f"Failed to fetch data. Status code: {response.status_code}")
    print('antibiotic data from STRING: ', len(lines))
#    print('anti data with needed columns: ', len(dict_list))
    return dict_list

def mapping(bacteria, antibiotics):
    """
    Maps PATRIC and STRING data together by taxon ID.
    bacteria, antibiotics: both list of dictionaries with taxon_id in their keys.
    """
    combined_data = []
    print("Combining PATRIC and STRING data into a combined list of dictionaries...")
    for bact in bacteria:
        match_found = False
        for anti in antibiotics:
            if bact['taxon_id'] == anti['#taxon_id']:
                combined_entry = {**bact, **anti}
                combined_data.append(combined_entry)
                match_found = True

    print('bacteria matched with STRING:', len(combined_data))
    return combined_data

def make_csv(data, file_name):
    """
    puts combined STRING and PATRIC data into a csv file in the current directory.
    data: list of dictionaries where each dictionary will get its own row.
    file_name: (str) name of csv file being created
    """
    existing_data = set()
    existing_fieldnames = None

    # Reading existing data
    try:
        with open(file_name, 'r', newline='') as input_file:
            dict_reader = csv.DictReader(input_file)
            existing_fieldnames = dict_reader.fieldnames
            for row in dict_reader:
                row_tuple = tuple(row.items())
                existing_data.add(row_tuple)
    except FileNotFoundError:
        print (f"{file_name} does not exist. a new file will be created.")

    new_data = []
    for row in data:
        # Standardize fieldnames
        standardized_row = {key.replace('#', ''): value for key, value in row.items()}
        row_tuple = tuple(standardized_row.items())
        if row_tuple not in existing_data:
            new_data.append(standardized_row)

    # If existing fieldnames aren't found, use fieldnames from new data
    if not existing_fieldnames and new_data:
        existing_fieldnames = new_data[0].keys()

    with open(file_name, 'a', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames=existing_fieldnames)

        if not existing_data:
            dict_writer.writeheader()

        dict_writer.writerows(new_data)

    print("Data has been written to", file_name)
#    print(f"Number of instances written to mapped_data.csv: {len(data)}")

def clear_csv(file_name):
    """
    Clears data written in a csv file.
    file_name: (str) name of csv
    """
    with open(file_name, 'w') as file:
        pass

def remove_incomplete_data(input_csv_path, column_name, output_csv_path):
    """
    Creates new csv after removing rows that are missing data in given column.
    
    input_csv_path: (str) path to csv path before removing rows
    column_name: (str) column where if the instance has no input in it, the row is removed.
    output_csv_path: (str) path to create the new csv file.
    """
    df = pd.read_csv(input_csv_path)
    if column_name not in df.columns:
        raise ValueError(f"Column {column_name} doesn't exist in input file.")
    cleaned_df = df.dropna(subset=[column_name])
    num_rows = len(cleaned_df)
    cleaned_df.to_csv(output_csv_path, index=False)

    print(f'Cleaned data saved to {output_csv_path}.')
#    print(f'Number of instances: {num_rows}')

def bacteria_data(data, url, folder_name='bacteria_proteins'):
    """
    Downloads bacteria protein data for graphs with needed taxon IDs and puts them in their
    own folder within the current directory.
    data: list of dictionaries with tax ID column
    url: (str) link to webpage with list of download links
    folder_name: (str) name of subfolder to create for downloads.
    """
    os.makedirs(folder_name, exist_ok=True)
    filter_values = {d['taxon_id'] for d in data}

    # Fetch page with download links
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch the page. Status code: {response.status_code}")
        return
    soup = BeautifulSoup(response.content, 'html.parser')
    links = soup.find_all('a', href=True)
    print(f"Fetched bacteria protein links from {url}")

    # Iterate over all links and matching taxon IDs
    all_links = 0
    mapped_links = 0
    skipped_links = 0
    for link in links:
        all_links += 1
        file_url = link['href']
        if any(file_url.startswith(value) for value in filter_values):
            mapped_links += 1
            if file_url.endswith(('.gz', '.size')):
                if not file_url.startswith('http'):
                    file_url = url.rstrip('/') + '/' +file_url.lstrip('/')
                    
                # Getting file names and preparing path
                filename = os.path.join(folder_name, os.path.basename(file_url))
                if filename.endswith('.gz'):
                    unzipped_filename = filename[:-3]
                elif filename.endswith('.size'):
                    unzipped_size_filename = filename[:-5]
                else:
                    unzipped_filename = filename
                    unzipped_size_filename = filename
                    
                # Check if the file already exists in folder
                if os.path.exists(filename) or os.path.exists(unzipped_filename) or os.path.exists(unzipped_size_filename):
                    skipped_links += 1
                    continue
                
                # Download the file
                file_response = requests.get(file_url, stream=True)
                if file_response.status_code == 200:
                    with open(filename, 'wb') as file:
                        for chunk in file_response.iter_content(chunk_size=8192):
                            file.write(chunk)
                    print(f"Downloaded {filename}")
                else:
                    print(f"Failed to download {file_url}. Status code: {file_response.status_code}")

#    print("links found for bacteria graphs", all_links)
#    print("links mapped from combined_data:", mapped_links)
#    print("Total links skipped (already existing):", skipped_links)
    print("All bacteria protein data downloaded.")

def sequences(data, url, folder_name='protein_sequences'):
    """
    Downloads protein sequences with taxon IDs that match the data and
    puts them in their own folder within the current directory.
    data: list of dictionaries with bacteria data with taxon ID column
    url: (str) link to webpage with list of download links
    folder_name: (str) name of subfolder to save downloaded protein data.
    """
    os.makedirs(folder_name, exist_ok=True)
    filter_values = {d['taxon_id'] for d in data}

    # Fetch page with download links
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch the page. Status code: {response.status_code}")
        return
    soup = BeautifulSoup(response.content, 'html.parser')
    links = soup.find_all('a', href=True)
    print(f"Fetched protein sequence data from {url}")

    # Iterate over all links
    all_links = 0
    mapped_links = 0
    skipped_links = 0
    for link in links:
        all_links += 1
        file_url = link['href']
        
        # Matching taxon IDs
        if any(file_url.startswith(value) for value in filter_values):
            mapped_links += 1
            if file_url.endswith(('.gz', '.size')):
                if not file_url.startswith('http'):
                    file_url = url.rstrip('/') + '/' +file_url.lstrip('/')
                    
                # Getting file names and preparing path
                filename = os.path.join(folder_name, os.path.basename(file_url))
                if filename.endswith('.gz'):
                    unzipped_filename = filename[:-3]
                elif filename.endswith('.size'):
                    unzipped_size_filename = filename[:-5]
                else:
                    unzipped_filename = filename
                    unzipped_size_filename = filename

                # Check if the file already exists in folder
                if os.path.exists(filename) or os.path.exists(unzipped_filename) or os.path.exists(unzipped_size_filename):
                    skipped_links += 1
                    continue
                
                # Downloading the file
                file_response = requests.get(file_url, stream=True)
                if file_response.status_code == 200:
                    with open(filename, 'wb') as file:
                        for chunk in file_response.iter_content(chunk_size=8192):
                            file.write(chunk)
                    print(f"Downloaded {filename}")
                else:
                    print(f"Failed to download {file_url}. Status code: {file_response.status_code}")

#    print("links found for protein sequences", all_links)
#    print("links mapped from combined_data:", mapped_links)
#    print("Total links skipped (already existing):", skipped_links)
    print("All protein sequence data downloaded.")

def unzip_files(folder_name):
    """
    Unzips all .gz files in specified folder
    folder_name: (str) name of folder containing .gz files
    """
    gz_files = []
    for root, _, files in os.walk(folder_name):
        print(f"Traversing directory: {root}")
        for file in files:
            if file.endswith('.gz') or file.endswith('.fa.gz'):
                gz_file_path = os.path.join(root, file)
                gz_files.append(gz_file_path)

    total_files = len(gz_files)
    if total_files == 0:
        print(f"No .gz files left in {folder_name}.")
        return

    print(f"Total .gz files to unzip in {folder_name}: {total_files}")

    # Process each .gz file
    for i, gz_file_path in enumerate(gz_files, 1):
        uncompressed_file_path = gz_file_path[:-3]
        if os.path.exists(uncompressed_file_path):
            print(f"{gz_file_path} is already unzipped.")
            continue
        try:
            with gzip.open(gz_file_path, 'rb') as f_in:
                with open(uncompressed_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(gz_file_path)
            print(f"Unzipped {gz_file_path}. Progress: {i}/{total_files}")
        except EOFError:
            print(f"EOFError: {gz_file_path} is corrupted and will be skipped.")
        except Exception as e:
            print(f"An error occurred while processing {gz_file_path}: {e}")

def unique_instances(csv_path, column_name):
    """
    Counts how many different values are in a column in a csv file.
    csv_file_path: (str) path to csv file to read
    column_name: (str) column name to check number of unique values
    """
    df = pd.read_csv(csv_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' does not exist in the CSV file.")
    unique_count = df[column_name].nunique()
    print(f'Number of unique {column_name}s:', unique_count)

def organize(parent_folder):
    """
    Puts all .gz.size files into a subfolder named 'gz_size_files'.
    parent_folder: (str) folder containing .gz.size files
    """
    subfolder_name = os.path.join(parent_folder, 'gz_size_files')
    os.makedirs(subfolder_name, exist_ok=True)
    print(f"Folder '{subfolder_name}' created or already exists.")

    for root, _, files in os.walk(parent_folder):
        for file in files:
            if file.endswith('.gz.size'):
                source_file = os.path.join(root, file)
                destination_file = os.path.join(subfolder_name, file)
                if os.path.exists(destination_file):
                    print(f"{file} already exists: skipping.")
                    continue
                shutil.move(source_file, destination_file)
                print(f"Moved {file} to {subfolder_name}")

    print("All .gz.size files have been moved to the subfolder in {parent_folder}")

def susceptibility(csv_path):
    """
    Counts and prints how many data instances are resistant vs. sustainable
    csv_path: (str) path to csv file containing data with 'resistant_phenotype' column
    """
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} doesn't exist for susceptibilitiy stats calculation.")
        return
    
    n_sus = 0
    n_res = 0
    column = 'resistant_phenotype'
    with open(csv_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            key = row.get(column, '').strip()
            if key == 'Susceptible':
                n_sus += 1
            elif key == 'Resistant':
                n_res += 1
    print("Number of susceptible bacteria:", n_sus)
    print ("Number of resistant bacteria:", n_res)

def all_graphs(input_folder, output_folder, delimiter=' ', max_workers=None):
    """
    Reads each bacteria protein download and calls make_graph to make a graph for the instance of data.
    All graphs are saved in one graph folder.
    input_folder: (str) path to folder containing .txt files
    output_folder: (str) path to the folder to save .gml graph files
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # creating a list taxon IDs of the graphs already created
    created_graphs = {filename.split('.')[0] for filename in os.listdir(output_folder) if filename.endswith('.gml')}
    count = len(created_graphs)
    print("created graphs:", count)

    # gathering list of files to process
    files_to_process = []
    for filename in os.listdir(input_folder):
        if filename.endswith('.txt'):
            taxon_id = filename.split('.')[0]
            if taxon_id not in created_graphs:
                files_to_process.append(os.path.join(input_folder, filename))
    to_do = len(files_to_process)
    print(f"{to_do} graphs to create...")

    # processing in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for filename in files_to_process:
            taxon_id = filename.split('.')[0]
            graph_path = os.path.join(output_folder, f"{taxon_id}.gml")
            file_path = os.path.join(input_folder, filename)
            futures.append(executor.submit(make_graph, file_path, graph_path, delimiter))
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error creating graph: {str(e)}")

    print("Done creating graphs!")

def make_graph(file_path, graph_path, delimiter):
    """
    Creates a graph from a single .txt file with nodes as the proteins and edges as their combined score.
    Graphs are created from their largest connected component.
    file_path (str): Path to the .txt file.
    graph_path (str): Path to save the .gml graph file.
    """
    G = nx.Graph()
    
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=delimiter)
        header = next(reader, None)
        for row in reader:
            if len(row) < 3:
                print(file_path, 'does not have required rows (skipping)')
                return
            else:
                protein1, protein2, combined_score = row[0], row[1], float(row[2])
                G.add_node(protein1)
                G.add_node(protein2)
                G.add_edge(protein1, protein2, weight=combined_score)

    largest = max(nx.connected_components(G), key=len)
    sub_graph = G.subgraph(largest).copy()
    nx.write_gml(sub_graph, graph_path)
    print(f"{graph_path} created.")

def statistics(folder_path, csv_path):
    """
    Creates a list of dictionaries with statistics on graphs in specified folder by calling process_graph,
    and saves the information to a csv file.
    folder_path: (str) folder name containing graphs saved in .gml format
    csv_path: (str) path to csv file to save statistics data
    """
        
    fieldnames = ['Taxon ID', 'Number of Nodes', 'Number of Edges',
            'Average Node Degree', 'Graph Density', 'Graph Diameter']

    with open(csv_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

    print("Gathering graph statistics...")
    stats = []
    file_paths = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if filename.endswith('.gml')]
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_graph, file_path): file_path for file_path in file_paths}
        
        for future in as_completed(futures):
            graph_stats = future.result()
            stats.append(graph_stats)

    graphs_processed = len(stats)
    print(f"Done with {graphs_processed} graphs.")

    with open(csv_path, mode='a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writerows(stats)
    print(f"statistics have been written to {csv_path}.")
    return stats

def process_graph(file_path):
    """
    Gathers statistics on a single .gml graph file.
    file_path: (str) path to .gml file
    """
    graph_stats = {}
    
    try:
        G = nx.read_gml(file_path)
        taxon_id = file_path[43:-4]
        
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        degree_dict = dict(nx.average_degree_connectivity(G))
        average_node_degree = sum(k * v for k, v in degree_dict.items()) / len(G)
        density = nx.density(G)
               
        if nx.is_connected(G):
            diameter = nx.diameter(G)
        else:
            diameter = "Graph is not connected"

        graph_stats['Taxon ID'] = taxon_id
        graph_stats['Number of Nodes'] = num_nodes
        graph_stats['Number of Edges'] = num_edges
        graph_stats['Average Node Degree'] = average_node_degree
        graph_stats['Graph Density'] = density
        graph_stats['Graph Diameter'] = diameter
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

    print(f"{file_path} stats done.")
    return graph_stats

def histograms(stats_path, output_folder):
    """
    Creates histograms from a csv file containing statistics on numerous graphs.
    stats_path: (str) path to csv with statistics data
    output_folder: (str) path to folder to save histograms showing stats
    """
    print("Creating histograms...")
    
    stats = []
    with open(stats_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            clean_row = {}
            for key, value in row.items():
                try:
                    float_value = float(value.strip())
                    clean_row[key] = float_value
                except ValueError:
                    clean_row[key] = None
            if any(value is not None for value in clean_row.values()):
                stats.append(clean_row)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    keys_to_plot = ['Number of Nodes', 'Number of Edges', 'Average Node Degree', 'Graph Density', 'Graph Diameter']

    for key in keys_to_plot:
        values = [entry[key] for entry in stats]
        
        plt.figure(figsize=(12, 9))
        plt.hist(values, bins=10, edgecolor='white', histtype='bar', color='skyblue')

        x_min, x_max = float(min(values)), float(max(values))
        plt.xlim(x_min * 0.8, x_max * 1.1)

        if key == 'Graph Diameter':
            plt.ylim = (0, 15)
        elif key == 'Average Node Degree':
            plt.ylim = (0, 300)
        else:
            plt.ylim = (0, 50)

        ax.set_xscale('log')
        plt.xlabel(f'{key}', fontsize=28)
        plt.ylabel('Frequency', fontsize=28)
        plt.title(f'{key} of Bacteria Graphs', fontsize=32)
        plt.xticks(fontsize=20)
        plt.tight_layout()

        output_file = os.path.join(output_folder, key)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"{key} histogram saved to {output_file}")

def graph_ratio(bacteria_proteins_path, subgraphs_path, folder_path):
    """
    Creates histograms comparing the size of the largest connected component
    graphs to their complete graphs.
    bacteria_proteins_path: (str) path to folder containing bacteria protein data
    used to create graphs
    subgraphs_path: (str) path to folder containing all create graphs from largest
    connect component
    folder_path: (str) path to folder to save histograms comparing size of graphs
    """
    # reading subgraphs csv
    subgraphs_list = []
    with open(subgraphs_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            subgraphs_list.append(row)
    
    needed_files = []
    for sub_graph in subgraphs_list:
        taxon_id = sub_graph.get('Taxon ID')
        needed_files.append(taxon_id)

    if not os.path.exists(bacteria_proteins_path):
        print('Main graphs path does not exist.')
        return

    total = []
    for filename in os.listdir(bacteria_proteins_path):
        taxon_id = filename.split('.')[0]
        if taxon_id in needed_files:
            single = {}
            file_path = os.path.join(bacteria_proteins_path, filename)

            # main graph
            nodes = set()
            main_edges = 0
            with open(file_path, 'r') as file:
                header = file.readline()
                for line in file:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split(' ')
                    if len(parts) < 3:
                        continue
                    node1, node2, _ = parts[:3]
                    nodes.add(node1)
                    nodes.add(node2)
                    main_edges += 1
            main_nodes = len(nodes)

            single['graph edges'] = main_edges
            single['graph nodes'] = main_nodes

            # subgraph
            taxon_id = filename.split('.')[0]
            for item in subgraphs_list:
                if item.get('Taxon ID') == taxon_id:
                    sub_edges = item.get('Number of Edges')
                    sub_nodes = item.get('Number of Nodes')
                    single['subgraph edges'] = sub_edges
                    single['subgraph nodes'] = sub_nodes
                    single['taxon_id'] = taxon_id
                    break
                    
            total.append(single)

    # ratio list of dictionaries
    ratios = []
    for graph in total:
        ratio = {}
        sub_nodes = int(graph.get('subgraph nodes'))
        sub_edges = int(graph.get('subgraph edges'))
        main_nodes = int(graph.get('graph nodes'))
        main_edges = int(graph.get('graph edges'))
        
        if main_nodes > 0:
            ratio_nodes = sub_nodes / main_nodes
        else:
            ratio_nodes = None

        if main_edges > 0:
            ratio_edges = sub_edges / main_edges
        else:
            ratio_edges = None

        ratio['Ratio of the Number of Nodes'] = ratio_nodes
        ratio['Ratio of the Number of Edges'] = ratio_edges
        ratios.append(ratio)

    # creating histograms
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    keys_to_plot = ['Ratio of the Number of Nodes', 'Ratio of the Number of Edges']

    for key in keys_to_plot:
        values = [entry[key] for entry in ratios]

        plt.figure(figsize=(12, 9))
        plt.hist(values, bins=10, edgecolor='white', histtype='bar', color='skyblue')
        if key == 'Ratio of the Number of Nodes':
            plt.xlabel('Ratio of Nodes', fontsize=28)
            plt.xlim(min(values), 1)
        else:
            plt.xlabel('Ratio of Edges', fontsize=28)
            plt.xlim(0.48, 0.51)
        plt.ylabel('Frequency', fontsize=28)
        plt.title(key, fontsize=32)
        plt.xticks(fontsize=20)
        plt.tight_layout()

        output_file = os.path.join(folder_path, key)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
def main():

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

if __name__ == "__main__":
    main()
