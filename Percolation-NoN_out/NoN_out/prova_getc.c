//
//  prova_getc.c
//  percolation_out
//
//  Created by gino on 2017-03-14.
//  Copyright © 2017 Gino. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>

#define MAX_LENGHT 1000000


int main(int argc, const char * argv[]) {

    int n;
    char c;
    char var[10000], *line, *start;
    
    line = malloc(MAX_LENGHT * sizeof(char));
    
    FILE *fp;
    fp= fopen("file.txt","r");
    
    
    while(fgets(line,MAX_LENGHT,fp) != NULL){
        start = line;
        printf("%s\n",line);
        while(sscanf(start,"%s%n",var,&n) == 1){
            printf("%s\n",var);
            start = start + n;
        }
    }
    
    
    
    
    
    return 0;
}
