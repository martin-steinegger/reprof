#ifndef UTIL_H
#define UTIL_H
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kseq.h"
#include "reprof_struct.h"
#include "blastpsimat.h"
KSEQ_INIT(gzFile, gzread)


#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
//                             A      C  D  E  F  G  H  I      K  L  M   N       P   Q   R   S   T       V   W   X   Y
const size_t lookup_aa[26] = { 0, -1, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, -1, 17, 18, 19, 20, -1 };   
const char revers_aa[26] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q','R', 'S', 'T', 'V', 'W', 'Y' };
size_t lookup_aa_to_index(char aa) {
    return lookup_aa[aa - 65];
}
// E  F  G   H  I   J   K  L
const size_t lookup_ss[8]  = { 0, -1, 1, 1, -1, -1, -1, 2};


size_t lookup_ss_to_index(char ss) {
    return lookup_ss[ss - 69];
}

const size_t NUMBER =0;
const size_t MASS=1;
const size_t VOLUME=2;
const size_t HYDROPHOBICITY =3;
const size_t CBETA = 4;
const size_t HBREAKER = 5; 
const size_t CHARGE = 6;
const size_t POLARITY = 7;

const char sec_converter[3] = {  'H', 'E', 'L' };



const float acc_norm[27] = {
    //A 
    106,  
    //B
    160,         
    //C
    135,  
    //D
    163, 
    //E 
    194,
    //F 
    197, 
    //G 
    84, 
    //H 
    184,
    //I 
    169, 
    //J 
    -1, 
    //K 
    205, 
    //L 
    164,
    //M 
    188, 
    //N 
    157,
    //O 
    -1, 
    //P 
    136,
    //Q 
    198, 
    //R 
    248, 
    //S 
    130,
    //T 
    142,
    //U 
    -1, 
    //V 
    142, 
    //W 
    227,
    //X 
    180,        
    //Y 
    222, 
    //Z 
    196,       
    //max
    248
};

const float aa_features[21][8] = { 
    //A C D E F G H I K L M N P Q R S T V W X Y
    //A 'number', 'mass', 'volume', 'hydrophobicity', 'cbeta', 'hbreaker', 'charge', 'polarity'
    { 0,   0.109, 0.170, 0.700,  0,   0,   0.5,     0},
    //C
    { 4,   0.357, 0.289,0.778, 0,    0,    0.5,      0},
    //D
    { 3,   0.450, 0.304,0.111, 0,    0,    0,        1},
    //E
    { 6,   0.558, 0.467,0.111, 0,    0,    0,        1},
    //F
    { 13,  0.698, 0.774,0.811, 0,    0,    0.5,      0},
    //G   
    { 7,   0.000, 0.000,0.456, 0,    0,    0.5,      0},
    //H
    { 8,   0.620, 0.555,0.144, 0,    0,    1,        1},
    //I
    { 9,   0.434, 0.636,1.000, 1,    0,    0.5,      0},
    //K
    { 11,  0.550, 0.647,0.067, 0,    0,    1,        1},
    //L
    { 10,  0.434, 0.636,0.922, 0,    0,    0.5,      0},
    //M
    { 12,  0.574, 0.613,0.711, 0,    0,    0.5,      0},
    //N
    { 2,   0.442, 0.322,0.111, 0,    0,    0.5,      1},
    //P
    { 14,  0.310, 0.314,0.322, 0,    1,    0.5,      0},
    //Q
    { 5,   0.550, 0.499,0.111, 0,    0,    0.5,      1},
    //R
    { 1,   0.767, 0.676,0.000, 0,    0,    1,        0},    
    //S
    { 15,  0.233, 0.172,0.411, 0,    0,    0.5,      1},
    //T
    { 16,  0.341, 0.334,0.422, 1,    0,    0.5,      1},
    //V
    { 19,  0.326, 0.476,0.967, 1,    0,    0.5,      0},
    //W
    { 17,  1.000, 1.000,0.400, 0,    0,    0.5,      0},
    //X
    { 20,  0,     0,    0,     0,    0,    0,        0},
    //Y
    { 18,  0.822, 0.796,0.356, 0,    0,    0.5,      1}
};


int number (const char aa) {
    return (int) aa_features[lookup_aa_to_index(aa)][NUMBER];
}

float mass (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][MASS];
}

float volume (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][VOLUME];
}

float hydrophobicity (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][HYDROPHOBICITY];
}
float cbeta (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][CBETA];
}
float hbreaker (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][HBREAKER];
}
float charge (const char aa) {
    return aa_features[lookup_aa_to_index(aa)][CHARGE];
}
float polarity (const char aa)  {
    return aa_features[lookup_aa_to_index(aa)][POLARITY];
}


float_array_2d * in_sequence_bit(size_t length) {
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** result;
    result = (float **) malloc(length*sizeof(float *)); 
    for(int i =0; i < length;i++){
        result[i] = (float *) malloc(1*sizeof(float));  
        result[i][0] = 1.0f;
    }
    ret_array->col_size=1;
    ret_array->row_size=length;
    ret_array->data = result;
    return ret_array;
}


char * parse_fasta(const char * fasta_file) {
    kseq_t *seq; 
    char * ret_seq;
    gzFile fp;
    fp = gzopen(fasta_file, "r");  
	seq = kseq_init(fp);
    if(kseq_read(seq) >= 0) {
        printf("name: %s\n", seq->name.s);  
        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);  
        printf("seq: %s\n", seq->seq.s); 
        ret_seq = malloc(strlen(seq->seq.s)*sizeof(char)+1);
        memcpy(ret_seq, seq->seq.s, strlen(seq->seq.s)*sizeof(char)+1);
        kseq_destroy(seq);
        gzclose(fp);
    }  
    return ret_seq;
}





float_array_2d * aa_composition (char * residues,size_t length) {
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** ret_composition;
    ret_composition = (float**) malloc(length*sizeof(float *));  
    float * composition;
    composition = (float*) calloc(21,sizeof(float));  
    
    for(int i =0; i <length;i++) {
        composition[number(residues[i])]++;
    }
    for(int i =0; i <21;i++) {
        composition[i] = composition[i] / length;
    }
    for(int i =0; i <length;i++) {
        ret_composition[i] = composition;
    }
    
    ret_array->col_size = 21;
    ret_array->row_size = length;
    ret_array->data = ret_composition;
    return ret_array;
}

float_array_2d * profile (char * residues,size_t length) {
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** profile;
    profile = (float**) malloc(length*sizeof(float*));  
    for (int i = 0; i < length; i++)  
        profile[i] = (float*) calloc(20,sizeof(float));  
    
    for(int i =0; i <length;i++) {
        int array_pos = number(residues[i]);
        if (array_pos < 20) {
            profile[i][array_pos] = 1.0f;
        }
    }
    ret_array->row_size = length;
    ret_array->col_size = 20;
    ret_array->data = profile;    
    return ret_array;
}

float_array_2d * distance_c(size_t length){
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** distances;
    
    distances = (float**) malloc(length*sizeof(float*));  
    for (int i = 0; i < length; i++)  
        distances[i] = (float*) malloc(4*sizeof(float));  
    for(int pre=1; pre <= length;pre++){
        size_t post = length - pre;
        for(int i=1; i <= 4; i++){
            if (post >= pow(2, (i - 1)) * 10) {
                distances[pre-1][i - 1] = 1;
            }
            else {
                float r = (post - pow(2,(i - 2)) * 10) / (pow(2, (i - 1)) * 10);
                
                distances[pre-1][i - 1] = (r < 0) ? 0: r;
            }
        }
    }
    ret_array->col_size = 4;
    ret_array->row_size = length;
    ret_array->data = distances;
    return ret_array;
}


float_array_2d * length_4state (size_t length) {
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** ret_fourstate=0;
    
    ret_fourstate = (float**) malloc(length*sizeof(float *));  
    float * result;
    result = (float*) calloc(0, (4*sizeof(float)));  
    
    for(int i=1; i <= 4; i++){
        if (length >= pow(2,(i - 1)) * 60) {
            result[i - 1] = 1;
        }
        else {
            float r = (length - pow(2, (i - 2)) * 60) / (pow(2,(i - 1)) * 60);
            result[i - 1] = (r < 0)? 0 : r;
        }
    }
    for(int i =0; i <length;i++) {
        ret_fourstate[i] = result;
    }
    ret_array->col_size = 4;
    ret_array->row_size = length;
    ret_array->data = ret_fourstate;
    return ret_array;
}


float_array_2d * distance_n(size_t length){
    float_array_2d * ret_array = malloc(sizeof(float_array_2d));
    float ** distances=0;
    
    distances = (float**) malloc(length*sizeof(float*));  
    for (int i = 0; i < length; i++)  
        distances[i] = (float*) malloc(4*sizeof(float));  
    for(int pre=0; pre < length;pre++){
        for(int i=1; i <= 4; i++){
            if (pre >= pow(2, (i - 1)) * 10) {
                distances[pre][i - 1] = 1;
            }
            else {
                float r = (pre - pow(2,(i - 2)) * 10) / (pow(2, (i - 1)) * 10);
                
                distances[pre][i - 1] = (r < 0) ? 0: r;
            }
        }
    }
    ret_array->col_size = 4;
    ret_array->row_size = length;
    ret_array->data = distances;
    return ret_array;
}



float_array_2d * calc_feature(char * feat,char * residues,size_t length){
    
    if (strcmp(feat, "distanceN") == 0) {
        static float_array_2d * distn = 0; 
        if(distn==0)
            distn=distance_n(length);
        return distn;
    } else if (strcmp(feat, "length_4state") == 0) {
        static float_array_2d * state = 0; 
        if(state==0)
            state=length_4state(length);
        return state;
    } else if (strcmp(feat, "profile") == 0) {
        static float_array_2d * prof = 0; 
        if(prof!=0){
            free_float_array_2d(prof);
            prof= profile(residues,length);
        }else{
            prof = profile(residues,length);
        }
        return prof;
    } else if (strcmp(feat, "aa_composition") == 0) {
        static float_array_2d * aa_comp = 0; 
        if(aa_comp!=0){
            free(aa_comp->data[0]);
            free_float_array_2d_ptr_only(aa_comp);
            aa_comp= aa_composition(residues,length);
        }else{
            aa_comp = aa_composition(residues,length);
        }
        return aa_comp;
    } else if (strcmp(feat, "distanceC") == 0) {
        static float_array_2d * distc = 0; 
        if(distc==0)
            distc=distance_c(length);
        return distc;
    } else if(strcmp(feat,"in_sequence_bit")==0){
        static float_array_2d * seq_bit = 0; 
        if(seq_bit==0)
            seq_bit=in_sequence_bit(length);
        return seq_bit;
    } else if(strcmp(feat,"normalized")==0){
        return &get_pis_blast_result()->normalized;
    } else if(strcmp(feat,"info")==0){
        return &get_pis_blast_result()->info;
    } else if(strcmp(feat,"weight")==0){
        return &get_pis_blast_result()->weight;
    } else if(strcmp(feat,"percentage")==0){
        return &get_pis_blast_result()->percentage;
    } else {
        return NULL;
    } 
}

size_t max_pos (float * array, size_t size){
    size_t mpos = 0;
    for(size_t pos =1; pos < size; pos++) {
        if (array[pos] > array[mpos]) {
            mpos = pos;
        }
    }
    return mpos;
}

float acc_rel2abs(float acc, char res) {
    return floor((acc / 100) * acc_norm[res-65]);
}

float acc_ten2rel(float * accs, size_t size) {
    size_t mpos = max_pos(accs,size);
    return floor( ((mpos * mpos) + ((mpos + 1) * (mpos + 1))) / 2 );
}

char acc_rel2three (float acc) {
    if (acc >= 36) 
        return 'e';
    else if (acc >= 9) 
        return 'i';
    else 
        return 'b';
}

char acc_rel2two (float acc) {
    if (acc >= 16) 
        return 'e';
    else 
        return 'b';
}

char sec_three2one(float * array, size_t size) {
    size_t mpos = max_pos(array,size);
    return sec_converter[mpos];
}

float reliability(float * array, size_t size) {
    size_t nr1 = max_pos(array,size);
    float val = array[nr1];
    float tmp=array[nr1];
    array[nr1] = 0;
    size_t nr2 = max_pos(array,size);
    float val2 = array[nr2];
    array[nr1]= tmp;
    return floor(10 * (val - val2));
}


float_array_2d * run_model(struct fann * model,float_array_list * inputs,size_t size){
    
    float_array_list current_input;
    float_array_2d * output;
    output = malloc(sizeof(float_array_2d ));
    output->data = (float **)malloc( size * sizeof(float *));
    output->col_size = model->num_output;
    output->row_size = size;
    for(size_t i = 0; i < size; i++){
        
        current_input=inputs[i];
        //       model->num_input = current_input.num_inuse;
        output->data[i] = (float *)malloc( model->num_output * sizeof(float ));
        memcpy(output->data[i],fann_run(model, current_input.list),model->num_output * sizeof(float ));
        free(inputs[i].list);
    }
    free(inputs);
    return output;
}

float_array_list * create_inputs( features_list * features, size_t from, size_t to
                                 ,char * sequ,size_t chain_length,
                                 size_t input_layer_size, float_array_2d * pre_output) {
    
    float_array_list * inputs=(float_array_list * )malloc((to-from+1)*sizeof(float_array_list));
    for(int i = 0; i < (to-from+1);i++){
        init_float_array_list(&inputs[i],input_layer_size);     
    }
    static float * float_null=0;
    if(float_null==0){
        float_null=(float*) calloc(chain_length,sizeof(float));  
    }
    float_array_2d * current_data;
    feature_item *feature=features;
    do{   // over all features 
        
        if (strcmp(feature->val->source, "output")==0) {
            
            current_data = pre_output;
        } else {
            if(strcmp(feature->val->source,"fasta")==0||
               strcmp(feature->val->source,"blastPsiMat")==0){ 
                current_data=calc_feature(feature->val->feature,sequ,chain_length);
            }else{
                continue;
            }
            //            current_data = $parsers{$source}->$feature;
        }
        
        int window=feature->val->window;
        for(int center = from; center <= to; center++){
            int win_start = center - (window - 1) / 2;
            int win_end  = center + (window - 1) / 2;   
            for(int iter=win_start;iter  <=win_end;iter++) {
                if (iter < 0 || iter >= chain_length) {
                    add_array_to_array(float_null, current_data->col_size, &inputs[center-from]);
                }
                else {    
                    add_array_to_array(current_data->data[iter], current_data->col_size, &inputs[center-from]);
                }
            }
        }
    }while((feature=feature->next)!=NULL);
    
    return inputs;
}   

void write_output(char * file_output, float_array_2d * sec_prediction, float_array_2d * acc_prediction,char * sequ,size_t chain_length){
    
    FILE *fp = fopen(file_output, "w");
    
    fprintf(fp,"##General\n");
    fprintf(fp,"# No\t: Residue number (beginning with 1)\n");
    fprintf(fp,"# AA\t: Amino acid\n");
    fprintf(fp,"##Secondary structure\n");
    fprintf(fp,"# PHEL\t: Secondary structure (H = Helix, E = Extended/Sheet, L = Loop)\n");
    fprintf(fp,"# RI_S\t: Reliability index (0 to 9 (most reliable))\n");
    fprintf(fp,"# pH\t: Probability helix (0 to 1)\n");
    fprintf(fp,"# pE\t: Probability extended (0 to 1)\n");
    fprintf(fp,"# pL\t: Probability loop (0 to 1)\n");
    fprintf(fp,"##Solvent accessibility\n");
    fprintf(fp,"# PACC\t: Absolute\n");
    fprintf(fp,"# PREL\t: Relative\n");
    fprintf(fp,"# P10\t: Relative in 10 states (0 - 9 (most exposed))\n");
    fprintf(fp,"# RI_A\t: Reliability index (0 to 9 (most reliable))\n");
    fprintf(fp,"# Pbe\t: Two states (b = buried, e = exposed)\n");
    fprintf(fp,"# Pbie\t: Three states (b = buried, i = intermediate, e = exposed)\n");
    fprintf(fp,"# \n");
    fprintf(fp,"No\tAA\tPHEL\tRI_S\tpH\tpE\tpL\tPACC\tPREL\tP10\tRI_A\tPbe\tPbie\n");
    for(size_t i=0; i < chain_length; i++) {
        size_t no = i + 1;
        char aa = sequ[i];
        
        char PHEL = sec_three2one(sec_prediction->data[i],sec_prediction->col_size);
        size_t RI_S = reliability(sec_prediction->data[i],sec_prediction->col_size);
        char PREL = acc_ten2rel(acc_prediction->data[i],acc_prediction->col_size);
        size_t PACC = acc_rel2abs(PREL, aa);
        char P10 = max_pos(acc_prediction->data[i],acc_prediction->col_size);
        char RI_A = reliability(acc_prediction->data[i],acc_prediction->col_size);
        size_t pH = floor(sec_prediction->data[i][0] * 100);
        size_t pE = floor(sec_prediction->data[i][1] * 100);
        size_t pL = floor(sec_prediction->data[i][2] * 100);
        char Pbe = acc_rel2two(PREL);
        char Pbie = acc_rel2three(PREL);
        fprintf(fp, "%u\t%c\t%c\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%c\t%c\n",
                no, aa, PHEL, RI_S, pH, pE, pL, PACC, PREL, P10, RI_A, Pbe, Pbie);
        //        say OUT join "\t", $No, $AA, $PHEL, $RI_S, $pH, $pE, $pL, $PACC, $PREL, $P10, $RI_A, $Pbe, $Pbie;
    }
    fclose(fp);
}

float sum_arr(float * myArray, size_t size)
{
    float sum = 0;
    for(size_t i = 0;i < size;i++)
    {
        sum = sum + myArray[i];
    }
    return(sum);
}


float_array_2d *  jury(float_array_2d * uu, float_array_2d *  ub, float_array_2d *  bu,float_array_2d *  bb){
    
    const float_array_2d * arrays[4] = { uu, ub, bu, bb}; 
    const size_t num_arrays = 4;
    const size_t num_pos = uu->col_size;
    //    
    
    float_array_2d * result=(float_array_2d *) malloc(1*sizeof(float_array_2d));
    
    result->col_size=num_pos; 
    result->row_size=uu->row_size;
    result->data = (float **) malloc(result->row_size*sizeof(float_array_2d));
    for(int i_line = 0; i_line < result->row_size; i_line++){
        float * tmp=(float*) calloc(num_pos,sizeof(float));  
        for(int i_array=0; i_array < num_arrays;i_array++) {
            float sum=sum_arr(arrays[i_array]->data[i_line],arrays[i_array]->col_size);
            for(int i_pos = 0; i_pos < num_pos ; i_pos++) {
                tmp[i_pos] += arrays[i_array]->data[i_line][i_pos] / sum;
                
            }
        }   
        for(int i = 0 ; i< num_pos;i++){
            tmp[i] = tmp[i] /num_arrays;
        }
        
        result->data[i_line]=tmp;
    }
    return result;
}



#endif