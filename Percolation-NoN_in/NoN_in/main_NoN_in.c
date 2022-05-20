//  main.c
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.

/* ***************************************************  */
/*   How to run: ./a.out N eps file_in_1 file_in_2 run  */
/* **************************************************** */

// where "run" is 1 if the code is run for the first time (first merging of the two clusters) or 2 if it's run after the first time (only one extra cluster has merged the Giant Component), N is the size of the network, eps is the threshold

#include "NoN_in.h"

#define directory_INPUT "/Users/ginodelferraro/Dropbox/GLASSO_Language/Glasso_Matlab/Jij_thresholded_e-4" //directory containing the input files


int main(int argc, const char * argv[]) {
    
    chdir(directory_INPUT); // Move to the Input directory
    
    struct variable *node_old, *node_new;       // pointer to a struct of type variable
    int *queue,**size_comp_old, **size_comp_new, **A_old, **A_new, tot_cl;  //  A adjacency matrix, tot_cl is tot numb of clusters
    double **J_old, **J_new;           // J input matrix
    
    input_parameters(argc,argv); // get size N, threshold eps
  
    read_allocate_Jij_Aij(argv,&J_old,&A_old,Input_1); // threshold Jij and create adjacency matrix Aij for Input file 1
    allocate_memory(&queue,&node_old); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
    count_size_clusters(node_old,&size_comp_old,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_old,2);
  
    
    read_allocate_Jij_Aij(argv,&J_new,&A_new,Input_2); // threshold Jij and create adjacency matrix Aij for Input file 2
    allocate_memory(&queue,&node_new); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm and label the nodes
    count_size_clusters(node_new,&size_comp_new,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_new,3);
    
    delta_size(node_old,node_new,size_comp_old,size_comp_new,A_old,run);
    
    free_memory(J_old,J_new,A_old,A_new,queue,size_comp_old,size_comp_new,node_old,node_new);

    
    return 0;
}





