//  main.c
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.

/* ***************************************************  */
/*   How to run: ./a.out N eps file_in_1 file_in_2 run  */
/* **************************************************** */

// where "run" is 1 if the code is run for the first time (first merging of the two clusters) or 2 if it's run after the first time (only one extra cluster has merged the Giant Component), N is the size of the network, eps is the threshold

#include "percolation_in.h"

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
    int *queue,**size_comp_old, **size_comp_new, **A_old, **A_new, **A_red, **A_middle;
    int **list_added, **list_tot, tot_cl, count,i,j,numb_cluster;
    double **J_old, **J_new, **J_middle;           // J input matrix
    
    count = 0;
    numb_cluster = 1;
    
    input_parameters(argc,argv); // get size N, threshold eps

    allocate_Ared_and_list(&A_red,&list_added,&list_tot);   // Allocate space for Adj. reduced, list of nodes added in module order, list of node added in nodes order
    
    run = 1;
    /* ---- BEFORE THE JUMP ----- */
 
    read_allocate_Jij_Aij(argv,&J_old,&A_old,Input_1); // threshold Jij and create adjacency matrix Aij for Input file 1, before the jump
    allocate_memory(&queue,&node_old);   // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
    count_size_of_clusters(node_old,&size_comp_old,tot_cl,N);  // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_old,node_old,GRAPH_1);
  
    
    /* ---- AFTER THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_new,&A_new,Input_2); // threshold Jij and create adjacency matrix Aij for Input file 2, after the jump
    allocate_memory(&queue,&node_new);   // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm and label the nodes
    count_size_of_clusters(node_new,&size_comp_new,tot_cl,N);   // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_new,node_old,GRAPH_2);
  //  free(A_new);
    
    load_A_frozen(node_old,node_new,size_comp_old,size_comp_new,A_old,A_red,list_added,list_tot,&count,run); // Add the freezed module(s) to the Adjacency matrix
    
    /* ---- AT THE MIDDLE OF THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_middle,&A_middle,Input_3); // threshold Jij and create adjacency matrix Aij for Input file 3, after the jump at the middle
    draw_graph(A_middle,node_old,GRAPH_3);
    

    compare_A_and_write(A_red,A_middle,list_added,list_tot,&count);
    
    printf("LT i |mod_i|\n");
    for (i=0; i < N ; i ++) {
        printf("%d %d\n",list_tot[i][0],list_tot[i][1]);
    }
    printf("LT i |mod_i|\n");
    for (i=0; i < count ; i ++) {
        printf("%d %d\n",list_added[i][0],list_added[i][1]);
    }
    
    free_memory(J_old,J_new,A_old,A_new,size_comp_old,size_comp_new,node_old,node_new);
    
    /* ------------------------------------ */
    
    
    run = 2;
    
    /* ---- BEFORE THE JUMP ----- */
    
    read_allocate_Jij_Aij(argv,&J_old,&A_old,Input_4); // threshold Jij and create adjacency matrix Aij for Input file 1, before the jump
    allocate_memory(&queue,&node_old);   // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
    count_size_of_clusters(node_old,&size_comp_old,tot_cl,N);  // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_old,node_old,GRAPH_4);
    
    
    /* ---- AFTER THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_new,&A_new,Input_5); // threshold Jij and create adjacency matrix Aij for Input file 2, after the jump
    allocate_memory(&queue,&node_new);   // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm and label the nodes
    count_size_of_clusters(node_new,&size_comp_new,tot_cl,N);   // Print Size of the 1st, 2nd and 3rd biggest cluster
    draw_graph(A_new,node_old,GRAPH_5);
    //  free(A_new);
    
    load_A_frozen(node_old,node_new,size_comp_old,size_comp_new,A_old,A_red,list_added,list_tot,&count,run); // Add the freezed module(s) to the Adjacency matrix
    
    /* ---- AT THE MIDDLE OF THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_middle,&A_middle,Input_6); // threshold Jij and create adjacency matrix Aij for Input file 3, after the jump at the middle
    draw_graph(A_middle,node_old,GRAPH_6);
    
    
    compare_A_and_write(A_red,A_middle,list_added,list_tot,&count);
    
    
    printf("LT i |mod_i|\n");
    for (i=0; i < N ; i ++) {
        printf("%d %d\n",list_tot[i][0],list_tot[i][1]);
    }
    printf("LT i |mod_i|\n");
    for (i=0; i < count ; i ++) {
        printf("%d %d\n",list_added[i][0],list_added[i][1]);
    }
   
   free_memory(J_old,J_new,A_old,A_new,size_comp_old,size_comp_new,node_old,node_new);
    
    
    return 0;
}





