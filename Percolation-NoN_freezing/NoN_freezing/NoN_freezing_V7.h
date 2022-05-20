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

#define MAX_LENGTH 1000


int N,num_jumps; // tot numb of nodes
double eps;

/* ------------------------*/
struct variable{  // structure for each node in the graph

    int LabCluster;    // Label of the node-cluster
    int SizeCluster;   // Size of the node-cluster
    int DeltaSize;   // Delta cluster Size between two lambdas
};

struct variable_vox{ // structure for each voxel coordinate and
    
    int node;
    int x,y,z; // x,y,z coordinates in scan
    double corr; // value of the correlation for that voxel

};


/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */


/*  get parameter from stdout and convert them in the right format */

void input_parameters(int argc, const char **argv, char ***file_jumps){
    
    int i;
    char filename[101];
    FILE *fp; // pointer to the file where the Jij names are stored
    
    if(argc != 4 ){
        printf("\nWrong number of input parameters! - EXIT \n\n");
        exit(errno);
    }

    if(argc == 4){
        
        N = atoi(argv[1]);    // Size of the network
        eps = atof(argv[2]);  // Threshold parameter
    }
    
    sprintf(filename,"%s",argv[3]); // file where the Jij names are stored
    fp = fopen(filename,"r");
    printf("file %s\n",filename);
    if(fp == NULL){
        fprintf(stderr,"PROBLEM OPENING INPUT FILE WITH Jij matrices names \n");
        exit(errno);
    }
    
    // create an array which contains the names of the files of the matrices Jij
    // for each jumps: before, after and at the plateau
    *file_jumps = (char **)malloc( 3 * sizeof(char *));
    for (i = 0; i < 3; i++) {
        (*file_jumps)[i] = (char *)malloc(101 * sizeof(char));
    }

    for(i = 0; i < 3; i++){
        
        fscanf(fp,"%s",(*file_jumps)[i]); // store the name of the matrices Jij into an array of char
        printf("%s\n",(*file_jumps)[i]);

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

void read_allocate_Jij_Aij(const char **argv, double ***J, int ***A, char **file_jumps, int file_in){
    
    int i,j,non_zeros;
    
   
    char filename[101];
    FILE *fp; // pointer to input file
    
    sprintf(filename,"%s",file_jumps[file_in]); // read name_file from stdin, file 1 or 2 depending on var
    fp=fopen(filename,"r");
    
    printf("Reading matrix %s\n",file_jumps[file_in]);
    
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
            
            if (fabs((*J)[i][j]) > eps){ // if |Jij| > eps, keep it and assign this value to Aij also
                
                (*A)[i][0] = (*A)[i][0] + 1;     // increment the degree of the node i, stored in A[i][0]
                (*A)[i][non_zeros] = j;     // Assign Jij to Aij
           
                non_zeros++; // increment numb of non zero elements

            }
            else (*J)[i][j] = 0.;  // if |Jij| < eps then threshold it to zero
            (*J)[i][i] = 0.; // set diagonal terms equal to zero
        }
    }
    
    return;

}

/* =========================================================================== */
/*                        FURTHER MEMORY ALLOCATION                            */
/* =========================================================================== */

void allocate_voxels(struct variable_vox **vox){
    
    *vox = (struct variable_vox *)calloc(N,sizeof(struct variable_vox)); // allocate space for the structure which contains x, y, z, corr of each node
    
    return;
}

/* =========================================================================== */

void allocate_memory(int **queue, struct variable **node){

    *queue = (int *)calloc(N,sizeof(int));  // Allocate memory for an array queue
    *node = (struct variable *)calloc(N,sizeof(struct variable));  // allocate memory for N struct variable
    
    return;
}

/* =========================================================================== */
/*                     Load voxel index, coordinates and correlations          */
/* =========================================================================== */


void load_voxels(struct variable_vox *vox){
    
    int i,module,n;
    char *line, *name, *start;
    
    line = calloc(1000, sizeof(char));
    name = calloc(1000, sizeof(char));
    
    FILE *fp;
    fp = fopen("sort_mod_relab.txt","r");
    
    if(fp == NULL){
        fprintf(stderr,"PROBLEM OPENING voxel file \n");
        exit(errno);
    }
    
    
    for (i = 0; i < N; i++) {
        
    
//          printf("%d %d %d %lf %d\n",vox[i].x,vox[i].y,vox[i].z,vox[i].corr,module);
        
        
        }
    
    
    
    printf("\n\n");
    
    fclose(fp);
    
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

    int i,j,temp,delta,neigh,Label_Cluster,label_1;
    int *current,*end;
    
    for (i=0; i< N; i++) {
        node[i].LabCluster = 0;
    }

    Label_Cluster = 0;
    
    label_1 = 0;
    
    for(i = 0; i < N; i++){
    
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
    
    printf("\ntot numb of clusters = %d\n",Label_Cluster);

    return Label_Cluster;  // Total number of clusters (how many)
}

/* ================================================================================ */
/*                         Count how many nodes per cluster                         */
/* ================================================================================ */

/* Count how many nodes there are in each cluster. Put in an array[M][2]            */
/* the label of the cluster and the size of it sorted descending based on the size  */


/* size_comp[i][0]is the label of the cluster */
/* size_comp[i][1]is the size of the cluster */
/* Note: labeling of the cluster start from 1 to M = tot_num_cl */

void count_size_of_clusters(struct variable *node, int ***size_comp, int tot_num_cl, int var, double *GC_before, int N){

    int i,j,LabClust,temp0,temp1;
    
    /* Allocate memory for the matrix which stores cluster-size and label */
    *size_comp = (int **)calloc((tot_num_cl+1),sizeof(int *));
    for (i=0; i<(tot_num_cl+1); i++)
        (*size_comp)[i] = (int *)calloc(2,sizeof(int));
    

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
    
    if(var == 0)
        *GC_before = (double)(*size_comp)[2][1]/N ; // Size of the GC before the jump
    
    if(tot_num_cl == 1){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
 //       printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t0\t0\t0\t0\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N);
    }
    
    if(tot_num_cl == 2){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
   //     printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t%d\t%lf\t0\t0\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N);
    }
    
    if(tot_num_cl == 3){
   //     printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N);
    }
    
    if(tot_num_cl == 4){
        //printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N,(*size_comp)[4][0],(double)(*size_comp)[4][1]/N);
    }
    
    if(tot_num_cl == 5){
     //   printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N,(*size_comp)[4][0],(double)(*size_comp)[4][1]/N,(*size_comp)[5][0],(double)(*size_comp)[5][1]/N);
    }
    
    if(tot_num_cl >= 6){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
        //  printf("%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N);
     //   printf("%c%c%c%c%c%c ",name_J_mat[p][2],name_J_mat[p][3],name_J_mat[p][4],name_J_mat[p][5],name_J_mat[p][6],name_J_mat[p][7]);
        printf("%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N,(*size_comp)[3][0],(double)(*size_comp)[3][1]/N,(*size_comp)[4][0],(double)(*size_comp)[4][1]/N,(*size_comp)[5][0],(double)(*size_comp)[5][1]/N,(*size_comp)[6][0],(double)(*size_comp)[6][1]/N);
    }
    
    if(var == 1)
        printf("--------- Size of the jump in the GC before and after the jump  = %lf\n",((double)(*size_comp)[2][1]/N - *GC_before));

    return ;

}





/* =========================================================================== */
/*                              Print nodes in modules                         */
/* =========================================================================== */

// This function writes the nodes and node clusters of the nodes belonging to the clusters which joined at the jump

void print_nodes(struct variable *node_old, struct variable *node_new, int **size_comp_new, const char **argv, struct variable_vox *vox){

    int i;
    
    char filename[101];
    FILE *fp, *fp_out;
    
    sprintf(filename,"%s",argv[4]);
    fp = fopen(filename,"w");
    
    fp_out = fopen("6th_jump.txt","w");
    

    for (i = 0; i < N; i++) {
        
        if((node_new[i].SizeCluster == size_comp_new[1][1]) && (node_old[i].SizeCluster > 5) ){ // if the label of the new node is equal to the label of GC and the cluster size is larger than a certan value
            
            fprintf(fp_out,"%d %d %d \n",i,node_old[i].LabCluster,node_old[i].SizeCluster);
            
            
        }
    }
    
    
    
    return;

}

/* =========================================================================== */
/*                              print 5 biggest clusters                       */
/* =========================================================================== */

void print_M_biggest_clusters(struct variable *node_old, int **size_comp_old, const char **argv){
    
    int i,cl, max_cl;
    max_cl = 5;
    
    char filename[101];
    FILE *fp;

    sprintf(filename,"%s",argv[5]);
    fp = fopen(filename,"w");
    
    
    for(i = 0; i < N; i++){
        
     //   printf("node %d\n",i);
        for(cl = 1; cl <= max_cl; cl++){
           // printf("label cl = %d, label comp = %d\n",node_old[i].LabCluster,size_comp_old[cl][0]);
            if(node_old[i].LabCluster == size_comp_old[cl][0]){
                fprintf(fp, "%d %d\n",i,node_old[i].LabCluster);
                break;
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
    
    return;
    
}

/* =========================================================================== */

void free_all(int **list_added,int **list_tot){

    free(list_added);
    free(list_tot);
    
    return;

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

















