#include "includes/fann.h"
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "util.h"
#include "reprof_struct.h"
#include "blastpsimat.h"
#include <sys/types.h>
#include <sys/stat.h>

//
//  main.c
//  madprof
// #--------------------------------------------------
// Predict secondary structure and solvent
// accessibility from sequece
//
// steinegger_martin@web.de
//  Created by Martin Steinegger on 25.03.12.
//  Copyright (c) 2012 -. All rights reserved.
//
features_list * parse_feature_file(const char *configuration_file);
void help();
int max_window=0;

const char *get_filename_ext(const char *filepath) {
    const char *dot = strrchr(filepath, '.');
    if(!dot || dot == filepath) return "";
    return dot + 1;
}

const char *get_filename(const char *filepath) {
    const char *filename = basename(filepath);
    const char *file = strtok(filename, ".");
    if(!file || file == filename) return filename;
    
    return file;
}

features_list * parse_feature_file(const char *configuration_file)
{
    feature_item * currElm, * head_list;
    head_list = NULL; 
	FILE *conf = fopen(configuration_file, "r");
	if(!conf)
	{
		return NULL;
	}    
    char line[200];
    
    while(fgets(line,sizeof(line),conf) != NULL){
        char head [80];
        struct feature * feat;
        feat = (struct feature *)malloc(sizeof(struct feature));
        sscanf(line,"%s ",&head); 
        if(strcmp(head,"option") == 0)
        {
            free(feat);
            continue;
            
        }
        if(strcmp(head,"output") == 0)
        {
            char windowStr [10];
            sscanf(line, "%s %s %s %s", &head,feat->source,feat->feature,windowStr);
            feat->window = atoi(windowStr);
            currElm = (feature_item *)malloc(sizeof(feature_item));
            currElm->val = feat;
            currElm->next  = head_list;
            head_list = currElm;
            continue;
        }
        if(strcmp(head,"input") == 0)
        {
            char windowStr [10];
            sscanf(line, "%s %s %s %s", &head,feat->source,feat->feature,windowStr);
            feat->window = atoi(windowStr);
            if (feat->window > max_window) {
                max_window = feat->window;
            }
            currElm = (feature_item *)malloc(sizeof(feature_item));
            currElm->val = feat;
            currElm->next  = head_list;
            head_list = currElm;
            continue;
        }
    }
	fclose(conf);
    feature_item *feature=head_list;
    currElm = NULL;
    feature_item *new_head_list=NULL;
    do{
        currElm = (feature_item *)malloc(sizeof(feature_item));
        currElm->val = (struct feature *)malloc(sizeof(struct feature));
        memcpy(currElm->val,feature->val,sizeof(struct feature));
        currElm->next  = new_head_list;
        new_head_list = currElm;
    } while((feature=feature->next)!=NULL);
    free_feature_list(head_list);
	return new_head_list;
}

void help(){
    printf("NAME:\n");
    printf("\tmadprof\n");
    printf("\n");
    printf("DESCRIPTION:\n");
    printf("    Secondary structure and solvent accessibility prediction using neural networsk\n");
    printf("    \n");
    printf("USAGE:\n");
    printf("    \tmadprof --fasta [query.fasta] --out [path to dir] --model [path to model]\n");
    printf("    \n");
    printf("OPTIONS:\n");
    printf("    \t--fasta\n");
    printf("    \tInput (single) FASTA file\n");
    printf("    \t--out\n");
    printf("    \tEither an output file or a directory. If not provided or a directory, the suffix of the input filename is replaced to create an output filename\n");
    printf("    \t--model\n");
    printf("    \tDirectory where the model and feature files are stored\n");
    
    printf("EXAMPLES:\n");
    printf("\n");
    printf("    \tmadprof -fasta query.fasta -out out_dir/\n");
}

int main (int argc, const char * argv[])
{
    
    const int model_size = 14;
    
    const int FA_MODEL = 0;
    const int FU_MODEL = 1;
    const int FB_MODEL = 2;
    const int FUU_MODEL = 3;
    const int FUB_MODEL = 4;
    const int FBU_MODEL = 5;
    const int FBB_MODEL = 6;
    const int A_MODEL = 7;
    const int U_MODEL = 8;
    const int B_MODEL = 9;
    const int UU_MODEL = 10;
    const int UB_MODEL = 11;
    const int BU_MODEL = 12;
    const int BB_MODEL = 13;
    
    const char model_names[14][4] = {"fa","fu", "fb", "fuu", "fub", "fbu", "fbb",
                                      "a", "u",  "b",  "uu",  "ub",  "bu",  "bb"};
    char pwd[256];
    getcwd(pwd,sizeof(pwd));
    size_t chain_length=0;
    char * sequ;
    char * out_file;
    const char * pathToModel="/usr/share/madprof/";
    const char * pathToOut=0;
    int hasPisMat=0;
    char outPahtBuf[256];
    int paramCount=0; 
    for (int i=1; i<argc; )
    {      
        
        if (!strcmp(argv[i],"--model")) {
            pathToModel=argv[++i];
        }else if (!strcmp(argv[i],"--out")||!strcmp(argv[i],"-o")){ 
            pathToOut=argv[++i];
			// FIXME: Error handling for unsuccessful case?
        } else if (!strcmp(argv[i],"--input")||!strcmp(argv[i],"-i")){
            
            if(!strcmp(get_filename_ext(argv[++i]),"blastPsiMat")){
                sequ=parse_blast_pis_mat(argv[i]);
                hasPisMat=1;
            } else {
                /* The MIT License start */
                sequ=parse_fasta(argv[i]);
                /* The MIT License end */
            }
            out_file=get_filename(argv[i]);
            paramCount++;
            chain_length = strlen(sequ);
        } else if(!strcmp(argv[i],"-h")) {help(); return 0;}
        else{
            i++;
        }
    } // end of for-loop for command line input
    if(argc == 1 )
    { help(); return 0;}
    if(paramCount!=1){
        printf("Error: Wrong parameter");
        return 0;
    }else{
        if(pathToOut==0){
            snprintf(outPahtBuf, sizeof outPahtBuf, "%s%s%s.reprof", pwd,"//",out_file);
        }else{
            snprintf(outPahtBuf, sizeof outPahtBuf, "%s", pathToOut);
        }
        
    }
    // Load models and feature lists
    struct fann **models;
    models = malloc(model_size * sizeof(struct fann* ));
    features_list **features;
    features = malloc(model_size * sizeof(struct features_in_out* ));
    
    for(int i  = 0; i<model_size; i++){
        char modelPahtBuf[256];
        char featurePahtBuf[256];
        snprintf(modelPahtBuf, sizeof modelPahtBuf, "%s%s%s%s", pathToModel,"//",(char*) &model_names[i],".model");
        snprintf(featurePahtBuf, sizeof featurePahtBuf, "%s%s%s%s", pathToModel,"//",(char*) &model_names[i],".features");
        models[i] = fann_create_from_file(modelPahtBuf);
        features[i] = parse_feature_file(featurePahtBuf);
    }
    
    
    
    //Precompute
    //Do prediction for the inputfile
    printf("seq main: %s %lu\n", sequ,chain_length); 
    if(hasPisMat==1){
        size_t input_size=(chain_length - 1) - 0+1;
        float_array_2d * a = run_model(models[A_MODEL], create_inputs(features[A_MODEL], 0, chain_length - 1,sequ,chain_length,models[A_MODEL]->num_input,NULL),input_size);
        float_array_2d * u = run_model(models[U_MODEL], create_inputs(features[U_MODEL], 0, chain_length - 1,sequ,chain_length,models[U_MODEL]->num_input,NULL),input_size);
        float_array_2d * b = run_model(models[B_MODEL], create_inputs(features[B_MODEL], 0, chain_length - 1,sequ,chain_length,models[B_MODEL]->num_input,NULL),input_size);
        
        float_array_2d * uu = run_model(models[UU_MODEL], create_inputs(features[UU_MODEL], 0, chain_length - 1,sequ,chain_length,models[UU_MODEL]->num_input, u),input_size);
        float_array_2d * ub = run_model(models[UB_MODEL], create_inputs(features[UB_MODEL], 0, chain_length - 1,sequ,chain_length,models[UB_MODEL]->num_input, u),input_size);
        float_array_2d * bu = run_model(models[BU_MODEL], create_inputs(features[BU_MODEL], 0, chain_length - 1,sequ,chain_length,models[BU_MODEL]->num_input, b),input_size);
        float_array_2d * bb = run_model(models[BB_MODEL], create_inputs(features[BB_MODEL], 0, chain_length - 1,sequ,chain_length,models[BB_MODEL]->num_input, b),input_size);
        
        float_array_2d * sec_ori=jury(uu, ub, bu, bb);
        
        write_output(&outPahtBuf, sec_ori, a, sequ, chain_length);
        free_float_array_2d(a);
        free_float_array_2d(u);
        free_float_array_2d(b);
        free_float_array_2d(uu);
        free_float_array_2d(ub);
        free_float_array_2d(bu);
        free_float_array_2d(bb);
        free_float_array_2d(sec_ori);
    }
    
    size_t input_size=(chain_length - 1) - 0+1;
    float_array_2d * a_ori = run_model(models[FA_MODEL], create_inputs(features[FA_MODEL], 0, chain_length - 1,sequ,chain_length,models[FA_MODEL]->num_input,NULL),input_size);
    float_array_2d * u_ori = run_model(models[FU_MODEL], create_inputs(features[FU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FU_MODEL]->num_input,NULL),input_size);
    float_array_2d * b_ori = run_model(models[FB_MODEL], create_inputs(features[FB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FB_MODEL]->num_input,NULL),input_size);
    
    float_array_2d * uu = run_model(models[FUU_MODEL], create_inputs(features[FUU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FUU_MODEL]->num_input, u_ori),input_size);
    float_array_2d * ub = run_model(models[FUB_MODEL], create_inputs(features[FUB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FUB_MODEL]->num_input, u_ori),input_size);
    float_array_2d * bu = run_model(models[FBU_MODEL], create_inputs(features[FBU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FBU_MODEL]->num_input, b_ori),input_size);
    float_array_2d * bb = run_model(models[FBB_MODEL], create_inputs(features[FBB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FBB_MODEL]->num_input, b_ori),input_size);
    
    float_array_2d * sec_ori=jury(uu, ub, bu, bb);
    char oriPahtBuf[256];
    snprintf(oriPahtBuf, sizeof oriPahtBuf, "%s_ORI", outPahtBuf);
    write_output(&oriPahtBuf, sec_ori, a_ori, sequ, chain_length);
    //Free memory
    free_float_array_2d(uu);
    free_float_array_2d(ub);
    free_float_array_2d(bu);
    free_float_array_2d(bb);
    
    
    for(int pos = 0; pos < chain_length; pos++){
        char aa_ori = sequ[pos];
        for(int mut = 0; mut < 20; mut++){
            if(aa_ori == revers_aa[mut])
                continue;
            char pathToSaveResult[256];
            sequ[pos] = revers_aa[mut];
            float_array_2d * sec = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
            memcpy(sec, sec_ori,1*sizeof(float_array_2d) );
            sec->data = malloc(sizeof (float*) * sec_ori->row_size);
            for(int i = 0; i < sec_ori->row_size; i++){
                sec->data[i] = malloc(sizeof (float) * sec_ori->row_size);
                memcpy(sec->data[i], sec_ori->data[i], sizeof (float) * sec_ori->col_size);
            }
            float_array_2d * a = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
            memcpy(a, a_ori,1*sizeof(float_array_2d) );
            a->data = malloc(sizeof (float*) * a_ori->row_size);
            for(int i = 0; i < a_ori->row_size; i++){
                a->data[i] = malloc(sizeof (float) * a_ori->row_size);
                memcpy(a->data[i], a_ori->data[i], sizeof (float) * a_ori->col_size);
            }
            
            float_array_2d * u = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
            memcpy(u, u_ori,1*sizeof(float_array_2d) );
            u->data = malloc(sizeof (float*) * u_ori->row_size);
            for(int i = 0; i < u_ori->row_size; i++){
                u->data[i] = malloc(sizeof (float) * u_ori->row_size);
                memcpy(u->data[i], u_ori->data[i], sizeof (float) * u_ori->col_size);
            }
            
            float_array_2d * b = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
            memcpy(b, b_ori,1*sizeof(float_array_2d) );
            b->data = malloc(sizeof (float*) * b_ori->row_size);
            for(int i = 0; i < b_ori->row_size; i++){
                b->data[i] = malloc(sizeof (float) * b_ori->row_size);
                memcpy(b->data[i], b_ori->data[i], sizeof (float) * b_ori->col_size);
            }
            
            
            size_t seq_from = max(((pos - 1) - 1 * max_window + 1), 0);
            size_t seq_to =   min(((pos - 1) + 1 * max_window - 1), (chain_length - 1));
            
            size_t struc_from = max(((pos - 1) - 2 * max_window + 2), 0);
            size_t struc_to   = min(((pos - 1) + 2 * max_window - 2), (chain_length - 1));
            
            
            input_size=seq_to- seq_from+1;
            float_array_2d * a_tmp = run_model(models[FA_MODEL], create_inputs(features[FA_MODEL], seq_from, seq_to,sequ,chain_length,models[FA_MODEL]->num_input,NULL),input_size);
            float_array_2d * u_tmp = run_model(models[FU_MODEL], create_inputs(features[FU_MODEL], seq_from, seq_to,sequ,chain_length,models[FU_MODEL]->num_input,NULL),input_size);
            float_array_2d * b_tmp = run_model(models[FB_MODEL], create_inputs(features[FB_MODEL], seq_from, seq_to,sequ,chain_length,models[FB_MODEL]->num_input,NULL),input_size);
            
            size_t iter = 0;
            for(size_t j=seq_from;j<= seq_to;j++,iter++) {
                free(a->data[j]);
                free(u->data[j]);
                free(b->data[j]);
                a->data[j] = a_tmp->data[iter];
                u->data[j] = u_tmp->data[iter];
                b->data[j] = b_tmp->data[iter];
            }
            
            input_size=struc_to- struc_from+1;
            float_array_2d * uu = run_model(models[FUU_MODEL], create_inputs(features[FUU_MODEL], struc_from, struc_to,sequ,chain_length,models[FUU_MODEL]->num_input, u),input_size);
            float_array_2d * ub = run_model(models[FUB_MODEL], create_inputs(features[FUB_MODEL], struc_from, struc_to,sequ,chain_length,models[FUB_MODEL]->num_input, u),input_size);
            float_array_2d * bu = run_model(models[FBU_MODEL], create_inputs(features[FBU_MODEL], struc_from, struc_to,sequ,chain_length,models[FBU_MODEL]->num_input, b),input_size);
            float_array_2d * bb = run_model(models[FBB_MODEL], create_inputs(features[FBB_MODEL], struc_from, struc_to,sequ,chain_length,models[FBB_MODEL]->num_input, b),input_size);
            float_array_2d * sec_tmp=jury(uu, ub, bu, bb);
            
            iter = 0;
            for(size_t j=struc_from;j<= struc_to;j++,iter++) {
                free(sec->data[j]);
                sec->data[j] = sec_tmp->data[iter];
            }
            //        
            
            printf( "%s%c%u%c\n", "Done: ",aa_ori,pos,revers_aa[mut]);
            snprintf(pathToSaveResult, sizeof pathToSaveResult, "%s%c%u%c", pathToOut,aa_ori,pos+1,revers_aa[mut]);
            
            write_output(pathToSaveResult, sec, a, sequ, chain_length);
            //free run
            free_float_array_2d(a);
            free_float_array_2d(u);
            free_float_array_2d(b);
            free_float_array_2d_ptr_only(a_tmp);
            free_float_array_2d_ptr_only(u_tmp);
            free_float_array_2d_ptr_only(b_tmp);
            
            free_float_array_2d(uu);
            free_float_array_2d(ub);
            free_float_array_2d(bu);
            free_float_array_2d(bb);
            free_float_array_2d(sec);
            free_float_array_2d_ptr_only(sec_tmp);
        }
        sequ[pos] = aa_ori;
    }
    
    
    
    for(int i = 0; i < sizeof(models) / sizeof(struct fann); i++){
        free(models[i]);
    }
    for(int i = 0; i < sizeof(features) / sizeof(features_list); i++){
        free_feature_list(features[i]);
    }
    return 0;
}





