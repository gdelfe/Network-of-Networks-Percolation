//  main.c
//  percolation_gino
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.

/* ***************************************************  */
/*   How to run: ./a.out N eps file_A_in file_A_tot         */
/* **************************************************** */


#include "NoN_out.h"

#define directory_INPUT "/Users/ginodelferraro/Dropbox/GLASSO_Language/Glasso_Matlab/Jij_thresholded_e-4" //directory containing the input files


int main(int argc, const char * argv[]) {
    
    chdir(directory_INPUT); // Move to the Input directory
    
    struct variable *node_old, *node_new;       // pointer to a struct of type variable
    int *queue, **size_comp_old, **size_comp_new, **A_in, **A_tot, lines;  //  A adjacency matrix, tot_cl is tot numb of clusters
    double **J_tot;           // J input matrix
    
    input_parameters(argc,argv); // get size N, threshold eps
  
    lines = read_allocate_Aij_in(argv,&A_in);
    
    read_allocate_Jij_Aij(argv,&J_tot,&A_tot); // threshold Jij and create adjacency matrix Aij for Input file 1

    
    compare_A_in_A_tot(A_in,lines,A_tot);
    
    /* allocate_memory(&queue,&node_old); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes

    count_size_clusters(node_old,&size_comp_old,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_old,1);
  

    
    read_allocate_Jij_Aij(argv,&J_new,&A_new,Input_2); // threshold Jij and create adjacency matrix Aij for Input file 2
    allocate_memory(&queue,&node_new); // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm and label the nodes
    count_size_clusters(node_new,&size_comp_new,tot_cl,N); // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_new,2);
    
    delta_size(node_old,node_new,size_comp_old,size_comp_new,A_old,run);
    
    free_memory(J_old,J_new,A_old,A_new,queue,size_comp_old,size_comp_new,node_old,node_new);
*/
    
    return 0;
}





