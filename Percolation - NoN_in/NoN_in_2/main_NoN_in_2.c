//
//  main.c
//  percolation_gino
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.
//

#include "NoN_in_2.h"

#define directory_INPUT "/Users/gino/Dropbox/GLASSO_Language/Glasso_Matlab/Jij_thresholded_e-4" //directory containing the input files


int main(int argc, const char * argv[]) {

    int i,j;
    
    chdir(directory_INPUT); // Move to the Input directory
    
    struct variable *node_old, *node_new;       // pointer to a struct of type variable
    int *queue,**size_comp, **A_old, **A_new, tot_cl;  //  A adjacency matrix, tot_cl is tot numb of clusters
    double **J_old,**J_new;           // J input matrix
    
    input_parameters(argc,argv); // get size N, threshold eps
  
    read_allocate_Jij_Aij(argv,&J_old,&A_old,Input_1); // threshold Jij and create adjacency matrix Aij for Input file 1
    allocate_memory(&queue,&node_old); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
    count_size_clusters(node_old,&size_comp,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_old,1);
  

    
    read_allocate_Jij_Aij(argv,&J_new,&A_new,Input_2); // threshold Jij and create adjacency matrix Aij for Input file 2
    allocate_memory(&queue,&node_new); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm and label the nodes
    count_size_clusters(node_new,&size_comp,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_new,2);
    
    delta_size(node_old,node_new,size_comp,A_old);
  
    
   // print_Jij_Aij();
   //
    // print_graph_component()
    
    free(J_old);
    free(J_new);
    free(A_old);
    free(A_new);
    free(queue);
    free(size_comp);
    free(node_old);
    free(node_new);
    
    
    return 0;
}
