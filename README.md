The Optimal Transport-based Network Null Model toolbox
is a series of MATLAB functions and scripts that generates and tests 
null models for complex networks.

The main directory is named OTNNM and contains the following sub-directories
-OT_FUNCTIONS that contains the main functions to generate null models and network sampling
-NETWORK_MODELS that contains functions to generate network with specific topological properties
-UTILITIES that contains auxiliary functions  
-PERFORMANCE_STUDIES that contains scripts for assessing the goodness of the OT approach in constraining the desired real network properties (in- and out-degree and strength and reciprocated links) while generating random network samples
-RESULTS_PROCESSING that contains scripts for generate plot of the results of the performance studies
-STUDIES_MAT_FILES that contains the output of the performance studies that are read by the scripts contained in the RESULTS_PROCESSING directory
-FIGURES that contains .png figures generated by scripts in the RESULTS_PROCESSING directory

OT_FUNCTIONS directory: 
 
[plan, totalCost] = sinkhorn_knopp(marginal1, marginal2, costMatrix, gamma, maxIter, tol)
Sinkhorn-Knopp Algorithm for Optimal Transport. This function computes the optimal transport plan with entropy regularization using the Sinkhorn-Knopp algorithm 

[plan, totalCost] = multi_obj_ot(marginal1, marginal2, costMatrix1, costMatrix2, gamma, lambda, max_iter)
Multi-Objective Optimal Transport. This function solves a multi-objective optimal transport problem using the Sinkhorn algorithm

[Ci, Q] = OTmodularity_dir(W, eta, plan)
Optimal Transport Modularity for Directed Networks. This function finds the optimal community structure and calculates modularity for a directed network

[X, Y, indices, marginal_X, marginal_Y, sampled_joint] = sample_2d_distribution(P, z)
Sample from a 2D Probability Distribution. This function samples points from a given 2D probability distribution
  

NETWORK_MODELS directory:

W = weighted_directed_graph(n, p, w_func)
Weighted Directed Graph Model. This function generates a weighted directed graph using a probabilistic approach

W = strength_driven_attachment(n, m, alpha)
Strength-Driven Attachment Model. This function generates a weighted undirected network using a strength-driven preferential attachment mechanism

W = fitness_based_model_directed(n, p, fitness_func, weight_func)
This function generates a directed weighted network based on node fitness

[W, communities] = community_directed_weighted_network(n, c, p_in, p_out, sparsity, asymmetry)
This function generates a directed weighted network with community structure
  

UTILITIES directory:

[in_strength, out_strength, in_degree, out_degree, rec_degree] =compute_network_properties(A)
Compute_network_properties. Calculates network properties from adjacency matrix

[Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options)
INPUTSDLG Enhanced input dialog box supporting multiple data types

[Ci,Q]=modularity_dir(W,eta)
Modularity for Directed Networks. This function finds the optimal community structure and calculates modularity for a directed network

[VIn] = variation_of_information(Cx, Cy)
VARIATION_OF_INFORMATION Calculates the normalized variation of information between community partitions

visualize_network(W, options, community_labels)
VISUALIZE_NETWORK Visualize a network given its adjacency matrix
    

PERFORMANCE_STUDIES directory:

WeightedRandomNetworkStudy.m
This script conducts a comprehensive study on weighted random directed networks.
It explores how various parameters affect network properties and compares the original network with null models generated using Optimal Transport (OT).

ReciprocatedStrengthDrivenAttachmentStudy.m
This script conducts a comprehensive study on reciprocated strength-driven attachment networks. It explores how various parameters affect network properties and compares the original network with null models generated using Multi-Objective Optimal Transport (OT).

StrengthDrivenAttachmentNetworkStudy.m
This script conducts a comprehensive study on strength-driven attachment networks.
It explores how various parameters affect network properties and compares the original network with null models generated using Optimal Transport (OT).
  
FitnessBasedModelDirectedStudy.m
This script conducts a comprehensive study on fitness-based directed network models.
It explores how various parameters affect network properties and compares the original network with null models generated using Optimal Transport (OT).
  
CommunityDirectedWeightedNetworksStudy.m
This script conducts a comprehensive study on community detection in directed weighted networks. It explores the performance of two community detection methods:
1) Standard modularity optimization
2) Optimal Transport (OT) based modularity optimization


RESULTS_PROCESSING directory:

ResultProcess.m
Network Model Analysis and Visualization Script. This script loads and analyzes results from different network models, allowing user selection of specific parameter combinations and visualizing the relationships between real and permuted network metrics.

ReciprocatedResultProcess.m
Network Analysis and Visualization Script. This script analyzes and visualizes the results of a network study focusing on reciprocated strength-driven attachment.

CommunityReciprocatedResultProcess.m
Network Parameter Analysis Script. This script analyzes and visualizes the effects of various parameters on network metrics in a complex network study. The metrics visualized include: 
1) Modularity for OT null model (OT) and Newman traditional modularity
2) Variation of information for OT null model (OT) and Newman traditional modularity   
# E-STATISTICAL_SOFTWARE-OTNNM
