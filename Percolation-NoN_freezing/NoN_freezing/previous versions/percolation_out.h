//  main.c
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

/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
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



/*******************************************************************************/
/* =========================================================================== */
/*                                    CORE                                     */
/* =========================================================================== */
/*******************************************************************************/

 /* Structure A_red = degree_i / j_1 / j_2 / etc ....  */

void compare_A_and_write(int **A_red, int **A_middle, int **list_added, int **list_tot, int *down, int *count){

    int i,j,k, nn, diff_module;
    
    FILE *fp;
    
    fp=fopen("Network.txt","w");
    
    if(fp==NULL){
        fprintf(stderr,"PROBLEM OPENING INPUT FILE Jij \n");
        exit(errno);
    }

    
    printf("i mod_i | deg_i | neighbours\n");
    for (i = 0; i < *count; i++){
        printf("%d %d | %d | ",list_added[i][0],list_added[i][1],A_red[i][0]);
        
        for (j = 1; j <= A_red[i][0]; j++){
            printf("%d ",A_red[i][j]);
        }
        printf("\n");
    }
        printf("\n");
    
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
    for(i = *down; i < *count; i++){
        
        k = list_added[i][0]; /* value of the node in the Adj reduce matrix */
        
        if(A_red[i][0] == A_middle[k][0]){ // If the degree is the same
            fprintf(fp,"deg same = %d %d    i = %d\n",A_red[i][0],A_middle[k][0],list_added[i][0]);

            for(j = 1; j <= A_red[i][0]; j++) // for all the neighbours
                fprintf(fp,"%d --> %d | %d %d | %d\n",list_added[i][0],A_red[i][j],list_added[i][1],list_added[i][1],0); // print i, j, mod_i, mod_j, type link = 1
        }
       
        if(A_red[i][0] != A_middle[k][0]){ // if the degree is different
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
    
    *down = *count;
 
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

















