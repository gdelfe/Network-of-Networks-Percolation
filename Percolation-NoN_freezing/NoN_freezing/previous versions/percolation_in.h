//
//  percolation_gino.h
//  percolation_gino
//
//  Created by Gino on 3/6/17.
//  Copyright (c) 2017 Gino. All rights reserved.
//

#ifndef percolation_gino_percolation_gino_h
#define percolation_gino_percolation_gino_h


#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>

#define RED red

int N,run; // tot numb of nodes
double eps;

/* ------------------------*/
struct variable{  // structure for each node in the graph

    int LabCluster;    // Label of the node-cluster
    int SizeCluster;   // Size of the node-cluster
    int DeltaSize;   // Delta cluster Size between two lambdas
};


/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */


/*  get parameter from stdout and convert them in the right format */

void input_parameters(int argc, const char **argv){
    
    
    if(argc != 9 ){
        printf("\nWrong number of input parameters! - EXIT \n\n");
        exit(errno);
    }

    if(argc == 9){
        
        N=atoi(argv[1]);    // Size of the network
        eps=atof(argv[2]);  // Threshold parameter
   //     run=atoi(argv[10]);  // See the headear of the main.c
    }

    return;
    
}

/* =========================================================================== */
/*                                  INPUT                                      */
/* =========================================================================== */

/* Read Jij from input file and store them into a matrix J[i][j], thresholding the elements < eps to zero. */
// Input file is passed as argument from stdin
// Create an adjacency matrix Aij where to store only the degree of i and the non-zero elements of Jij
// Structure of Aij: A[i] =  degree_i, element_1, element_2 etc..

void read_allocate_Jij_Aij(const char **argv, double ***J, int ***A, int file_in){
    
    int i,j,non_zeros;
    
   
    char filename[101];
    FILE *fp; // pointer to input file
    
    sprintf(filename,"%s",argv[file_in]); // read name_file from stdin, file 1 or 2 depending on var
    fp=fopen(filename,"r");
    
    if(fp==NULL){
        fprintf(stderr,"PROBLEM OPENING INPUT FILE\n");
        exit(errno);
    }
    
    
    /* ALLOCATE MEMORY FOR  J[i][j] and A[i][j] ========================= */
    
    *J = (double **)malloc( N * sizeof(double *));
    *A = (int **)malloc( N * sizeof(int *));
    
    for( i = 0; i < N; i++){
        (*J)[i] = (double *)malloc( N * sizeof(double));
        (*A)[i] = (int *)malloc(N * sizeof(int)); // Allocate only one element per row, i.e. the space dedicated to store the degree of node i
        (*A)[i][0] = 0.; // set the degree of node i to 0
        
    }
    /* ==============================================================  */
    
    
    // read and threshold Jij -------- assign non-zero values to Aij
    for (i = 0; i < N; i++){
        
        non_zeros = 1;
        for (j=0; j < N; j++){ // read one line
         fscanf(fp, "%lf",& (*J)[i][j]); // read element
            
            if ((*J)[i][j] > eps){ // if Jij is larger then eps, keep it and assign this value to Aij also
                
                (*A)[i][0] = (*A)[i][0] + 1;     // increment the degree of the node i, stored in A[i][0]
                (*A)[i][non_zeros] = j;     // Assign Jij to Aij
           
                non_zeros++; // increment numb of non zero elements

            }
            else (*J)[i][j] = 0.;  // if Jij < eps the threshold it to zero
        }
    }
    
    
    return;

}

/* =========================================================================== */
/*                        FURTHER MEMORY ALLOCATION                            */
/* =========================================================================== */

void allocate_memory(int **queue, struct variable **node){

    *queue = (int *)calloc(N,sizeof(int));  // Allocate memory for an array queue
    *node = (struct variable *)calloc(N,sizeof(struct variable));  // allocate memory for N struct variable
    
    return;
}

/************************************************************************************/
/* =========================================================================== */
/*                                    CORE                                     */
/* =========================================================================== */
/************************************************************************************/



/* =========================================================================== */
/*                         BREADTH FIRST SEARCH                                */
/* =========================================================================== */

/* Breadth First Search algorithm. Label each node with the cluster label of belonging */

int compute_graph_components(struct variable *node, int **A, int *queue, int N){

    int i,j,temp,delta,neigh,Label_Cluster;
    int *current,*end;
    
    for (i=0; i< N; i++) {
        node[i].LabCluster = 0;
    }

    Label_Cluster = 0;
    
    for(i=0; i < N; i++){
    
        if(node[i].LabCluster == 0){ // if the node has not been already visited
        
            Label_Cluster++ ;      // increment the size of the cluster
            
            node[i].LabCluster = Label_Cluster;
            
            queue[0] = i;
            
            current = queue;    // pointer to the current element in the queue
            end = queue + 1;    // pointer to the last element added in the queue
            
            temp = 1;
            delta = 1;
            
            while (current != end) {  // while queue is not empty
                for (j=1; j<= A[*current][0]; j++) { // for all the neighbourds of i
                    
                    neigh = A[*current][j];
                    
                    if (node[neigh].LabCluster == 0) {
                        
                        node[neigh].LabCluster = Label_Cluster;
                        
                        queue[temp] = neigh; //assign the neighbour to the queue
                        temp++;
                        
                    }
                }
                
                current += 1; // point to next node in the queue
                end += temp - delta; // point to the end of the queue
                delta = temp;
            }
        }
    }
    free(queue);
    
    printf("tot clusters = %d\n",Label_Cluster);

    return Label_Cluster;  // Total number of clusters (how many)
}

/* ================================================================================ */
/*                         Count how many nodes per cluster                         */
/* ================================================================================ */

/* Count how many nodes there are in each cluster. Put in an array[M][2]            */
/* the label of the cluster and the size of it sorted descending based on the size  */
 
void count_size_of_clusters(struct variable *node, int ***size_comp, int tot_num_cl, int N){

    int i,j,LabClust,temp0,temp1;
    
    /* Allocate memory for a mat[M][2] where M is the total number of different clusters */
    /* size_comp[i][0]is the label of the cluster */
    /* size_comp[i][1]is the size of the cluster */
    /* Note: labeling of the cluster start from 1 to M = tot_num_cl */
    
    *size_comp = (int **)calloc((tot_num_cl+1),sizeof(int *));
    for (i=0; i<(tot_num_cl+1); i++)
        (*size_comp)[i] = (int *)calloc(2,sizeof(int));
    
    for (i=0; i< (tot_num_cl+1); i++) {
        (*size_comp)[i][0] = 0;  // set label cluster to 0
        (*size_comp)[i][1] = 0;  // set size cluster to 0
    }

    /* Count the size of each cluster */
    for (i = 0; i < N; i++){

        LabClust = node[i].LabCluster;
        
        (*size_comp)[LabClust][0] = node[i].LabCluster;  // Label of the cluster
        (*size_comp)[LabClust][1] = (*size_comp)[LabClust][1] + 1;  // count how many node are in each cluster
        }
    
    
     /* Assign the size of each cluster to the variable in the structure */
    for (i = 0; i < N; i++){
        node[i].SizeCluster = (*size_comp)[node[i].LabCluster][1];
    }
    
    /* sort the vector size_comp in descending order based on the size of the clusters, second column  */
    /* BUBBLE SORT */
    for (i = 1; i < tot_num_cl+1; i++){
        for (j = i+1; j < tot_num_cl+1;j++){
            if ((*size_comp)[i][1] < (*size_comp)[j][1])
            {
                temp1 = (*size_comp)[i][1];
                (*size_comp)[i][1] = (*size_comp)[j][1];
                (*size_comp)[j][1] = temp1;
                
                temp0 = (*size_comp)[i][0];
                (*size_comp)[i][0] = (*size_comp)[j][0];
                (*size_comp)[j][0] = temp0;
            }
        }
    }

    if(tot_num_cl==2){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
        printf("%d\t%lf\t%d\t%lf\t0\t0\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N);
    }
    
    if(tot_num_cl==3){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
       printf("%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N);
    }
    /* Print Label GC, Size GC */
   // else  printf("%d\t%lf\t0\t0\t0\t0\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N);

    return ;

}


/* =========================================================================== */
/*                              PRINT GRAPH                                    */
/* =========================================================================== */

/* Print a graph in gephViz format. Read the graph and plot it in a pdf file */
/* then remove the data file containing the graph */

void draw_graph(int **A, struct variable *node_old, int var){

    int i,j;
    
    FILE *fp;
    char filename[101];
    char color[5][100];
    
    sprintf(color[0],"gold");
    sprintf(color[1],"cyan");
    sprintf(color[2],"brown1");
    sprintf(color[3],"grey");
    sprintf(color[4],"ghostwhite");
    
    
    sprintf(filename,"graph_%d.dat",var); // read name_file from stdin, file 1 or 2 depending on var
    fp = fopen(filename,"w");
    
  
    fprintf(fp,"strict graph {\n");
    fprintf(fp,"    node[style=filled]\n");
    for (i=0; i<N; i++) { //print different colors for different clusters
        fprintf(fp,"    node[fillcolor=\"%s\"] %d\n",color[node_old[i].LabCluster-1],i);
    }
     for (i=0; i<N; i++) { // print edges
       if(A[i][0]==0) fprintf(fp," %d\n",i);
        for (j=1; j<= A[i][0]; j++) {
            fprintf(fp," %d -- %d\n",i,A[i][j]);
        }
    }
    fprintf(fp,"}\n");
    fclose(fp);
    
    sprintf(filename,"dot -Tpdf graph_%d.dat -o graph_%d.pdf",var,var);
    system(filename);
    sprintf(filename,"rm graph_%d.dat",var);
    system(filename);

    return ;

}




/* =========================================================================== */
/*               Construct the Adjacency matrix made of frozen clusters        */
/* =========================================================================== */

/* Load the frozen clusters and the INLINK connections as they were BEFORE the jump */

void load_A_frozen(struct variable *node_old, struct variable *node_new, int **size_comp_old, int **size_comp_new, int **A, int **A_red, int **list_added, int **list_tot, int *count,int run){

    int i,j;
    
    FILE *fp;
//    char filename[101];
  //  sprintf(filename, "cluster",var);
    
    if(run==1) system("rm Adj_freezed.txt");
    
    fp = fopen("Adj_freezed.txt","a");
    
    if(fp == NULL){
        fprintf(stderr,"PROBLEM OPENING OUTPUT FILE\n");
        exit(errno);
    }
    
    
    // print nodes, clusters new and old
    for (i=0; i<N; i++) {
        node_new[i].DeltaSize = node_new[i].SizeCluster - node_old[i].SizeCluster;
        printf("i = %d      Cl_old = %d   S_old = %d      Cl_new = %d   S_new = %d      diff Size = %d\n",i,node_old[i].LabCluster,node_old[i].SizeCluster,node_new[i].LabCluster,node_new[i].SizeCluster,node_new[i].DeltaSize);
    }
    printf("GC label = %d GC size = %d \n",size_comp_new[1][0],size_comp_new[1][1]);
    
   
    /* If the Label of the node is equal to the Label of the GC, then it has been merged with GC */
    /* Add it to the Adjacency matrix the node, cluster and neighbours as it was before merging */
    /* Structure of output A_in = cluster_i / degree_i / i / j_1 / j_2 / etc ....         */
    
    fprintf(fp,"i mod_i | deg_i | neighbours \n");
    
    if(run == 1){    /* If the first two clusters are merging for the first time into the GC put both in the Adjacency matrix */
        for (i = 0; i < N; i++) {
            
            if(node_new[i].LabCluster == size_comp_new[1][0]){ // if the label of the new node is equal to the label of GC
                
                A_red[*count][0] = A[i][0]; /* degree i */
                list_added[*count][0] = i;  /* add node i to the added_node_list consecutively */
                list_added[*count][1] = node_old[i].LabCluster; /* add label of i */
                list_tot[i][0] = i; /* add node i to the list in position i */
                list_tot[i][1] = node_old[i].LabCluster;    /* add label of i in position i */
                
                fprintf(fp, "%d %d | %d | ",i,node_old[i].LabCluster,A[i][0]);  // Print to file the Adj matrix of the nodes which joined the GC: cluser_i / deg_i / i / j_1, j_2, etc...
                for (j=1; j <= A[i][0]; j++) {
                    
                    A_red[*count][j] = A[i][j];
                    fprintf(fp,"%d ",A[i][j]);
                }
                fprintf(fp,"\n");
                (*count)++;       // increment the number of nodes added to A_reduced
            }
        }
    }
    
    if(run == 2){  // If the first two clusters already merged, add to NoN_in only the new cluster which joined

       for (i = 0; i<N; i++) {
            if(node_new[i].LabCluster == size_comp_new[1][0] && node_old[i].LabCluster != size_comp_old[1][0] ){ // If the node is now in the GC but it wasn't before
               
                A_red[*count][0] = A[i][0]; // degree i
                list_added[*count][0] = i;   //add node i to the added_node_list consecutively
                list_added[*count][1] = node_old[i].LabCluster; // add label of i
                list_tot[i][0] = i; // add node i to the list in position i
                list_tot[i][1] = node_old[i].LabCluster;    // add label of i in position i
                
                fprintf(fp, "%d %d | %d | ",i,node_old[i].LabCluster,A[i][0]);  // Print to file the Adj matrix of the nodes which joined the GC: cluser_i / deg_i / i / j_1, j_2, etc...
                for (j=1; j <= A[i][0]; j++) {
                    
                    A_red[*count][j] = A[i][j];
                    fprintf(fp,"%d ",A[i][j]);
                }
                fprintf(fp,"\n");
                (*count)++;       // increment the number of nodes added
            
            }
        }
    }
    
    
    fclose(fp);

    return ;
}


/* =========================================================================== */
/*           Allocate memory for A reduced and list of nodes added             */
/* =========================================================================== */


void allocate_Ared_and_list(int ***A_red,int ***list_added,int ***list_tot){
    
    int i;
    
    *A_red = (int **)malloc(N * sizeof(int *));
    *list_added = (int **)malloc(N * sizeof(int *));
    *list_tot = (int **)malloc(N * sizeof(int *));
    
    for (i = 0; i < N; i++) {
        
        (*A_red)[i] = (int *)malloc(N * sizeof(int));
        (*list_added)[i] = (int *)malloc(2 * sizeof(int));
        (*list_tot)[i] = (int *)calloc(2,sizeof(int));
        (*list_tot)[i][1] = -1;
    }
    
    return;
}


/* =========================================================================== */
/*       Compare A frozen with A at the middle and write OUTLINKS              */
/* =========================================================================== */


/* Structure A_red = degree_i / j_1 / j_2 / etc ....  */

void compare_A_and_write(int **A_red, int **A_middle, int **list_added, int **list_tot, int *count){
    
    int i,j,k, nn, diff_module;
    
    FILE *fp;
    
    fp=fopen("Network.txt","w");
    
    if(fp==NULL){
        fprintf(stderr,"PROBLEM OPENING INPUT FILE Jij \n");
        exit(errno);
    }
    
    // Matrix A_red contains only the frozen clusters and the connections among them BEFORE the jump
    printf("i mod_i | deg_i | neighbours\n");
    for (i = 0; i < *count; i++){
        printf("%d %d | %d | ",list_added[i][0],list_added[i][1],A_red[i][0]); // print: i  mod_i | degree_i
        
        for (j = 1; j <= A_red[i][0]; j++){ // print the neighbours
            printf("%d ",A_red[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    // Matrix A_middle contains all the connections AFTER the jump, at the middle of the "plateau"
    printf("i mod_i | deg_i | neighbours\n");
    for (i = 0; i < N; i++){
        printf("%d %d | %d | ",list_tot[i][0],list_tot[i][1],A_middle[i][0]);
        
        for (j = 1; j <= A_middle[i][0]; j++){
            printf("%d ",A_middle[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    printf("\n");
    
    
    /****************/
    
    fprintf(fp,"i j | mod_i mod_j | TL \n");
    for(i = 0; i < *count; i++){
        
        k = list_added[i][0]; /* value of the node in the Adj reduced matrix */
        
        if(A_red[i][0] == A_middle[k][0]){ // If the degree is the same -- just write the neighbours
            fprintf(fp,"deg same = %d %d    i = %d\n",A_red[i][0],A_middle[k][0],list_added[i][0]);
    
            for(j = 1; j <= A_red[i][0]; j++) // for all the neighbours
                fprintf(fp,"%d --> %d | %d %d | %d\n",list_added[i][0],A_red[i][j],list_added[i][1],list_added[i][1],0); // print i, j, mod_i, mod_j, type link = 1
        }
        
        // if the degree is different -- check if the new nodes are part of the frozen cluster added
        // or if they have joined after the jump, in this latter case do not add them
        if(A_red[i][0] != A_middle[k][0]){
            fprintf(fp,"deg diff = %d %d    i = %d\n",A_red[i][0],A_middle[k][0],list_added[i][0]);
            
            for(j = 1; j <= A_middle[i][0]; j++){ // for all the neighbours
                
                nn = A_middle[k][j];  // pick the neighbour in the Adj matrix middle
                // printf("neigh = %d\n",nn);
                if(list_tot[nn][1]!= -1){ // if the neighbours is in one of the previous clusters -- add it
                    
                    diff_module = (list_tot[k][1])!=(list_tot[nn][1]); // check if the two nodes are in the same cluster. Result: 1 if different, 0 if the same
                    
                    fprintf(fp,"%d --> %d | %d %d | %d\n",list_tot[k][0],list_tot[nn][0],list_tot[k][1],list_tot[nn][1],diff_module); // print i, j, mod_i, mod_j, type link = 1
                }
            }
        }
        
    }
    
    
    return;
    
}

/* =========================================================================== */
/*                              FREE ALLOCATED MEMORY                          */
/* =========================================================================== */


void free_memory(double **J_old,double **J_new,int **A_old,int **A_new,int **size_comp_old,int **size_comp_new,struct variable *node_old, struct variable *node_new){
    
    free(J_old);
    free(J_new);
    free(A_old);
    free(A_new);
    
    free(size_comp_old);
    free(size_comp_new);
    free(node_old);
    free(node_new);
    
    
}



































/************************************************************************************/
/* QUICKSORT function for a matrix Nx2 where the sorting is made based on the first column */

/** Divide  : Partition the array A[low....high] into two sub-arrays
 *           A[low....j-1] and A[j+1...high] such that each element
 *           of A[low....j-1] is less than or equal to A[j], which
 *           in turn is is less than or equal to A[j+1...high]. Compute
 *           the index j as part of this partitioning procedure.
 * Conquer : Sort the two sub-arrays A[low....j-1] and A[j+1....high]
 *           by recursive calls to quicksort
 **/
/*
void quick_sort(int **mat, int low, int high){
    
    int pivot, i, j, temp0, temp1;
    if(low < high) {
        pivot = low; // select a pivot element
        i = low;
        j = high;
        
        while(i < j) {
            // increment i till you get a number greater than the pivot element
            while(mat[i][0] <= mat[pivot][0] && i <= high)
                i++;
            // decrement j till you get a number less than the pivot element
            while(mat[j][0] > mat[pivot][0] && j >= low)
                j--;
            // if i < j swap the elements in locations i and j
            if(i < j) {
                temp0 = mat[i][0];
            //    temp1 = mat[i][1];
                
                mat[i][0] = mat[j][0];
                mat[j][0] = temp0;
                
            //    mat[i][1] = mat[j][1];
            //    mat[j][1] = temp1;
                
            }
        }
        
        // when i >= j it means the j-th position is the correct position
        // of the pivot element, hence swap the pivot element with the
        // element in the j-th position
        temp0 = mat[j][0];
        mat[j][0] = mat[pivot][0];
        mat[pivot][0] = temp0;
        
     //   temp1 = mat[j][1];
      //  mat[j][1] = mat[pivot][1];
      //  mat[pivot][1] = temp1;
        
        // Repeat quicksort for the two sub-arrays, one to the left of j
        // and one to the right of j
        quick_sort(mat, low, j-1);
        quick_sort(mat, j+1, high);
    }
}


*/
/*
void quick_sort(int **mat, int low, int high){

    int left,right,pivot,temp;
    
    if (high>1){
    
        pivot = low;
        left = low;
        right = high;
        
        while (left <= right) {
            while (mat[left][1] < pivot)
                left++;
            while (mat[right][1] > pivot)
                right--;
            if (left <= right) {
                
                temp = mat[right][1];
                mat[right][1] = mat[left][1];
                mat[left][1] = temp;
                
                left++;
                right--;
            }
        }
        quick_sort(mat,0,right);
        quick_sort(mat,left,high);
    
    }
    
    return;

}

*/

















