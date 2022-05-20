//  main.c
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.

/* ***************************************************  */
/*   How to run: ./a.out N eps file_in_1 file_in_2 jump  */
/* **************************************************** */

// where "jumo" is 1 if the code is run for the first time (first merging of the two clusters) or 2 if it's run after the first time (only one extra cluster has merged the Giant Component), N is the size of the network, eps is the threshold

#include "percolation_in_v3.h"

#define Input_1 3
#define Input_2 4
#define Input_3 5
#define Input_4 6
#define Input_5 7
#define Input_6 8

#define GRAPH_1 1
#define GRAPH_2 2
#define GRAPH_3 3
#define GRAPH_4 4
#define GRAPH_5 5
#define GRAPH_6 6

#define directory_INPUT "/Users/ginodelferraro/Dropbox/GLASSO_Language/Glasso_Matlab/Jij_thresholded_e-4/output" //directory containing the input files


int main(int argc, const char * argv[]) {
    
    chdir(directory_INPUT); // Move to the Input directory
    
    struct variable *node_old, *node_new;       // pointer to a struct of type variable
    int *queue,**size_comp_old, **size_comp_new, **list_added, **list_tot, Input_file, Graph_out,i,j;
    int **A_old, **A_new, **A_red, **A_middle,**A_middle_tot, tot_cl, count, count_cl;
    double **J_old, **J_new, **J_middle;           // J input matrix
    
    
    count = 0;          // tot number of node added in the final network
    count_cl = 3;      // number of clusters added in the final network
    Input_file = 2;
    Graph_out = 0;
    
    input_parameters(argc,argv); // get size N, threshold eps

    allocate_A_tot_and_list(&A_red,&A_middle_tot,&list_added,&list_tot);   // Allocate space for Adj. reduced, list of nodes added in module order,
                                                                              // list of node added in nodes order
    
    for (jump = 1; jump <= 2 ; jump++){
    
        /* ---- BEFORE THE JUMP ----- */
        read_allocate_Jij_Aij(argv,&J_old,&A_old,++Input_file); // threshold Jij and create adjacency matrix Aij for Input file 1, before the jump
        allocate_memory(&queue,&node_old);                      // allocate memory for queue and structure node_old[i]
        tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
        count_size_of_clusters(node_old,&size_comp_old,tot_cl,N);  // Print Size of the 1st, 2nd and 3rd biggest cluster
        
        
        
        /* ---- AFTER THE JUMP ----- */
        read_allocate_Jij_Aij(argv,&J_new,&A_new,++Input_file); // threshold Jij and create adjacency matrix Aij for Input file 2, after the jump
        allocate_memory(&queue,&node_new);                      // allocate memory for queue and structure node_new[i]
        tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm
        count_size_of_clusters(node_new,&size_comp_new,tot_cl,N);   // Print Size of the 1st, 2nd and 3rd biggest cluster
        
    
        
        /* ---- AT THE MIDDLE OF THE PLATEAU ----- */
        read_allocate_Jij_Aij(argv,&J_middle,&A_middle,++Input_file); // threshold Jij and create Aij for Input file 3, at the plateau after the jump
        
        
        /* ---- CREATE A-FROZEN AND A-MIDDLE-TOT ---- */
        /* Add the freezed module(s) before the jump to the Adjacency matrix A_red and at the plateau to Adj_middle_tot */
        /* Re-label the clusters in 1,2,3 order of addition */
        
        load_A_tot(node_old,node_new,size_comp_old,size_comp_new,A_old,A_red,A_middle,A_middle_tot,list_added,list_tot,&count,&count_cl,jump);
        //  free(A_new);
        
    /*    draw_graph(A_old,node_old,++Graph_out); // before the jump
        draw_graph(A_new,node_old,++Graph_out);  // right after the jump
        draw_graph(A_middle,node_old,++Graph_out); // at the plateau
      */  
        
        free_memory(J_old,J_new,A_old,A_new,size_comp_old,size_comp_new,node_old,node_new);

    }
    
    /* Compare A red with A middle tot and write the final Adj matrix */
    compare_A_and_write(A_red,A_middle_tot,list_added,list_tot,&count);
    
    /* draw the resulting final graph */
    draw_graph_final(&count,list_added);

    free_all(A_red,A_middle_tot,list_added,list_tot); // free memory
    
    
    return 0;
}






/*
 
 printf("i LT |mod_i|\n");
 for (i=0; i < N ; i ++) {
 printf("%d %d\n",list_tot[i][0],list_tot[i][1]);
 }
 printf("i LT |mod_i|\n");
 for (i=0; i < count ; i ++) {
 printf("%d %d\n",list_added[i][0],list_added[i][1]);
 }
 
 */


