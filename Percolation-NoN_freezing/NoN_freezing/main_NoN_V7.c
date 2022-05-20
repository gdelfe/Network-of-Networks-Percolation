//  main.c
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.

/* ***************************************************  */
/*   How to run: ./a.out N eps file_J_matrices file_out_nodes */ /** NOTE: file_J_matrices should be in the input folder (see below) **/
/* **************************************************** */

// where "jump" is 1 if the code is run for the first time (first merging of the two clusters) or 2 if it's run after the first time (only one extra cluster has merged the Giant Component), N is the size of the network, eps is the threshold

#include "NoN_freezing_V7.h"

#define directory_INPUT "/Volumes/LAB/BRAIN/DATA/MSKCC_raw_data/Healthy/4/J_C_data" //directory containing the input files
#define directory_OUTPUT "/Volumes/LAB/BRAIN/DATA/MSKCC_raw_data/Healthy/4/NoN" //directory for the output files

#define BEFORE 0
#define AFTER 1 

int main(int argc, const char * argv[]) {
    
    struct variable *node_old, *node_new;       // pointer to a struct of type variable
    struct variable_vox *vox;               // pointer to a struct of type voxel
    int *queue,**size_comp_old, **size_comp_new,i,j;
    int **A_old, **A_new, tot_cl;
    double **J_old, **J_new, GC_before;           // J input matrix
    char **file_jumps;

    
    chdir(directory_INPUT); // Move to the Input directory
    
    input_parameters(argc,argv,&file_jumps); // get size N, threshold eps, matrices at the jump
    allocate_voxels(&vox);
    
    load_voxels(vox);
   
    /* ---- BEFORE THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_old,&A_old,file_jumps,BEFORE); // threshold Jij and create adjacency matrix Aij for Input file 1, before the jump
    allocate_memory(&queue,&node_old);                      // allocate memory for queue and structure node_old[i]
    tot_cl = compute_graph_components(node_old,A_old,queue,N);  // Run BFS algorithm and label the nodes
    count_size_of_clusters(node_old,&size_comp_old,tot_cl,BEFORE,&GC_before,N);  // Print Size of the 1st, 2nd and 3rd biggest cluster
    
    
    /* ---- AFTER THE JUMP ----- */
    read_allocate_Jij_Aij(argv,&J_new,&A_new,file_jumps,AFTER); // threshold Jij and create adjacency matrix Aij for Input file 2, after the jump
    allocate_memory(&queue,&node_new);                      // allocate memory for queue and structure node_new[i]
    tot_cl = compute_graph_components(node_new,A_new,queue,N);  // Run BFS algorithm
    count_size_of_clusters(node_new,&size_comp_new,tot_cl,AFTER,&GC_before,N);   // Print Size of the 1st, 2nd and 3rd biggest cluster
    
    
    chdir(directory_OUTPUT); // Move to the OUTPUT directory
   
    print_nodes(node_old,node_new,size_comp_new,argv,vox);
    
   // print_M_biggest_clusters(node_old,size_comp_old,argv);
    
    
    free_memory(J_old,J_new,A_old,A_new,size_comp_old,size_comp_new,node_old,node_new);
    
    printf("Writing output\n");

    
    
    return 0;
}






/*
 
 printf("i LT |mod_i|\n");
 for (i=0; i < N ; i ++) {
 printf("%d %d\n",list_tot[i][0],list_tot[i][1]);
 }
c
 
 */


