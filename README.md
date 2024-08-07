# Leveraging Graph Convolutional Networks to Predict Antibiotic Susceptibility to Bacterial Strains.

  This research aims to develop a predictive model using graph convolutional networks (GCNs) to accurately predict the effectiveness of antibiotics against unique bacterial strains. The GCN model will be trained on an extensive dataset comprised of thousands of experimental microbiology measurements. By leveraging this data, the model will learn the intricate patterns and relationships between bacterial genetics and their antibiotic susceptibility. This can significantly enhance the understanding of bacterial resistance mechanisms and fost the development of new antibiotics and treatment strategies. Further, it can aid the identification of suitable antibiotics for the treatment of bacterial infections, thus improving clinical outcomes and combatting the growth of antibiotic resistance.

Node Features and Drug Association:
  Each instance of bacteria data is represented by its own graph, with nodes denoting proteins, and edges symbolizing protein-protein interactions. Each node will have features generated from amino acid sequences of the bacterial proteins, and some will have drug association scores that indicate how strongly an antibiotic targets a protein. However, it's unlikely all bacterial proteins will have that information for multiple antibiotics.

Data Collection and Curation:
  Large scale data on antimicrobial resistance (AMR) was obtained from the PATRIC database. From this database, instances comprising of bacterial strain, antibiotic, and susceptibility information were extracted. Then, these bacteria strains were mapped to the STRING database of protein-protein interaction networks. Node embeddings were then calculated from amino acid sequence data also from STRING. Targets for antibiotics and their corresponding associated scores will be extracted from the STITCH database of interaction networks of chemicals and proteins.

GCN Model Structure:
- 3 convolutional layers that will use the ReLU activation function to capture hierarchical structures in protein interaction networks.
- A global pooling layer that will distill node features into a graph representation
- 2 fully connected layers that will use the ReLU activation function to refine the model's features.
- An output layer that will use the Softmax function to predict whether a strain is sensitive or resistant to an antibiotic
- Binary cross-entropy loss function and Adam optimizer will be used for training.
- Five-fold cross validation, evaluating effectiveness, accuracy, precision, recall, F1-score, and ROC-AUC score.

Repository Contents:
  This algorithm written in Python, generates the graphs to later be used in training the GCN model. It includes the process of data curation, obtaining the AMR data, the protein sequence data, and the protein-protein interaction data. It also assesses the diversity of the curated data and graphs generated from that data, producing histograms using Networkx assessing the number of nodes, number of edges, graph density, graph diameter, and average degree of nodes in each created graph. The number of unique bacterial strains and the number of resistant vs. susceptible data is also found.

Acknowledgements:
This material is based upon work supported by the National Science Foundation under award OAC-2150491 with additional support from the Center for Computation & Technology at Louisiana State University, and my research mentor Michal Brylinski.
Computer support is provided by HPC@LSU computing.
