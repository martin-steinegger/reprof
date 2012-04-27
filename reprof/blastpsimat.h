//
//  blastpsimat.c
//  reprof
//
//  Created by Martin Steinegger on 29.03.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#ifndef BLASTPSIMAT_H
#define BLASTPSIMAT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "reprof_struct.h"


pis_blast_mat * get_pis_blast_result(){
    static pis_blast_mat * ret_val=0;
    if(ret_val==0)
        ret_val = malloc(sizeof(pis_blast_mat));
    return ret_val;
}

float normalize(float x) {
    return 1.0 / (1.0 + exp(-x));
}


char * parse_blast_pis_mat(const char *blastPsiMat) {
    FILE *conf = fopen(blastPsiMat, "r");
	if(!conf)
	{
		return NULL;
	}    
    int c=0,b;while ((b=fgetc(conf))!=EOF) c+=(b==10)?1:0;fseek(conf,0,SEEK_SET);
    pis_blast_mat * result = get_pis_blast_result();
    result->sequ=malloc(sizeof(char)*c);
    result->raw.data=malloc(sizeof(float *) *c);
    result->raw.col_size=20;
    result->percentage.data=malloc(sizeof(float *) *c);
    result->percentage.col_size=20;
    result->normalized.data=malloc(sizeof(float *) *c);
    result->normalized.col_size=20;
    result->info.data=malloc(sizeof(float *) *c);
    result->info.col_size=1;
    result->weight.data=malloc(sizeof(float *) *c);
    result->weight.col_size=1;
    char line[200];
    int data=0;
    int currPos=0;
    while(fgets(line,sizeof(line),conf) != NULL){
        char * pch;
        pch = strtok (line," ");
        if(!strcmp(pch,"Last")){
            data=1;
            continue;
        }
        if(data==1){
            if(isdigit(pch[0])){
                pch = strtok (NULL, " ");
                result->sequ[currPos]=pch[0];
                result->raw.data[currPos]=malloc(sizeof(float)*20);
                result->normalized.data[currPos]=malloc(sizeof(float)*20);
                for(int i =0; i < 20;i++){
                    pch = strtok (NULL, " ");
                    result->raw.data[currPos][i]=atoi(pch);
                    result->normalized.data[currPos][i]=normalize(result->raw.data[currPos][i]);
                }
                result->percentage.data[currPos]=malloc(sizeof(float)*20);
                for(int i =0; i < 20; i++){
                    pch = strtok (NULL, " ");
                    float numb = atoi(pch);
                    result->percentage.data[currPos][i]=numb/100.0f;
                }
                pch = strtok (NULL, " ");
                result->info.data[currPos]=malloc(sizeof(float)*1);
                result->info.data[currPos][0]=atof(pch);
                pch = strtok (NULL, " ");
                result->weight.data[currPos]=malloc(sizeof(float)*1);
                if (strncmp(pch,"inf",3)==0) {
                    result->weight.data[currPos][0]=0.0f;
                }
                else {
                    result->weight.data[currPos][0]=atof(pch);
                }


                currPos++;
                continue;
            }
        }
    }
    size_t chain_length = strlen(result->sequ);
    result->raw.row_size=chain_length;
    result->percentage.row_size=chain_length;
    result->normalized.row_size=chain_length;
    result->info.row_size=chain_length;
    result->weight.row_size=chain_length;
    return result->sequ;
}

#endif
