//
//  reprof_struct.c
//  reprof
// #--------------------------------------------------
// Predict secondary structure and solvent
// accessibility from sequece
//
// steinegger_martin@web.de
//  Created by Martin Steinegger on 25.03.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#ifndef REPROF_STRUCT_H
#define REPROF_STRUCT_H
#include <stdio.h>




struct feature {
    char source[80];
    char feature[80];
    int  window;
};
struct struct_features_list {
    struct feature * val;
    struct struct_features_list * next;
};
typedef struct struct_features_list feature_item;
typedef struct struct_features_list features_list;

struct struct_float_array_list
{
    size_t    num_alloc;
    size_t    num_inuse;
    float   *list;
};
typedef struct struct_float_array_list float_array_list;

struct struct_float_array_2d {
    size_t row_size;
    size_t col_size;
    float ** data;
};
typedef struct struct_float_array_2d float_array_2d;


struct struct_pis_blast_mat {
    char * sequ;
    float_array_2d raw;
    float_array_2d normalized;
    float_array_2d percentage;
    float_array_2d info;
    float_array_2d weight;
};
typedef struct struct_pis_blast_mat pis_blast_mat;

void print_float_array_list(float_array_list * toprint,size_t sizeof_toprint);
void free_feature_list(features_list *head);
void init_float_array_list(float_array_list *array,size_t init_size);
void add_array_to_array(float * list_to_add, size_t size,float_array_list * array);
void free_float_array_list(float_array_list * array);
void free_float_array_2d(float_array_2d * array);
void free_float_array_2d_ptr_only(float_array_2d * array);
void print_float_array_2d(float_array_2d * array);

void print_float_array_list(float_array_list * toprint,size_t sizeof_toprint){
    printf("[\n");
    for(int col = 0; col < sizeof_toprint; col++){
        printf("  [\n");
        for(int i = 0 ; i < toprint->num_inuse; i++){
            printf("        %f, \n", toprint[col].list[i]);
        }
        printf("  ]\n");
    }
    printf("]\n");
}


void free_feature_list(features_list *head){
    feature_item *tmp;
    while (head != NULL) {
        free(head->val);    
        tmp = head->next;
        free(head);
        head = tmp;
    }
}




void init_float_array_list(float_array_list *array,size_t init_size)
{
    array->num_inuse = 0;
    array->num_alloc = init_size;
    array->list      = malloc(init_size*sizeof(float));
}

void add_array_to_array(float * list_to_add, size_t size,float_array_list * array)
{
    if ((array->num_inuse+size) > array->num_alloc)
    {
        size_t new_size = ((array->num_alloc + size) * 2);
        float *new_list = realloc(array->list, new_size * sizeof(float));
        if (new_list == 0)
            printf("out of memory");
        array->num_alloc = new_size;
        array->list      = new_list;
    }
    memcpy(array->list+(array->num_inuse), list_to_add, sizeof(float)*size);
    array->num_inuse =  array->num_inuse+size;
}

void free_float_array_list(float_array_list * array){
    free(array->list);
    free(array);
}


void free_float_array_2d_ptr_only(float_array_2d * array){
    free(array->data);
    free(array);
}

void free_float_array_2d(float_array_2d * array){
    for(int row=0; row < array->row_size;row++)
        if(array->data[row]!=0) { 
            free(array->data[row]);
        }
    free(array->data);
    free(array);
}

void print_float_array_2d(float_array_2d * array){
    printf("[\n");
    for(int row=0; row< array->row_size;row++){
        printf("  [\n");        
        for(int col=0; col < array->col_size; col++){
            printf("       %f,\n", array->data[row][col]);
        }
        printf("  ]\n");
    }
    printf("]\n");
    
}
#endif
