#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim_bp.h"
#include <math.h>
/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/
#define bit_extraction(addr,n,k) (addr >> (n-k))
#define calc_index(addr,m) (((1 << (m)) - 1) & (addr >> 2)) 
#define calc_history_reg(addr,m) (((1 << (m)) - 1) & (addr)) 
#define calc_size(index) (int)(pow(2,index))

void initialise_predictor_table(predictor_table_t *this,unsigned long int index_bimodal,int k,unsigned long int index_gshare){
	int i;
	this->table_index = index_bimodal;
	this->bk = (block_t *) malloc ((calc_size(index_bimodal))*(sizeof(block_t)));
	for (i=0; i < calc_size(index_bimodal); i++){
		this->bk[i].counter = 2;
	}
    this->bk_h = (block_t *) malloc ((calc_size(index_gshare))*(sizeof(block_t)));
	for (i=0; i < calc_size(index_gshare); i++){
		this->bk_h[i].counter = 2;
	}
    this->bk_ct = (block_t *) malloc ((calc_size(k))*(sizeof(block_t)));
	for (i=0; i < calc_size(k); i++){
		this->bk_ct[i].counter = 1;
	}

	this->global_history_register = 0;
	this->predictions_values.miss_prediction = 0;
	this->predictions_values.prediction_event = 0;
}

void prediction(predictor_table_t *this,unsigned long int index,char outcome){
	this->predictions_values.prediction_event++;
	if((outcome == 't' && this->bk[index].counter < 2) | (outcome == 'n' && this->bk[index].counter > 1)) this->predictions_values.miss_prediction++;

	if(outcome == 't' && this->bk[index].counter != 3) this->bk[index].counter++;
	else if(outcome == 'n' && this->bk[index].counter != 0) this->bk[index].counter--;
}

void gshare_prediction(predictor_table_t *this,unsigned long int index,int m,int n, char outcome){
	unsigned long int n_part = bit_extraction(index,m,n);
	unsigned long int n_part_history_register = this->global_history_register;
	
	unsigned long int gshare_index = n_part_history_register ^ n_part;
	unsigned long int final_index = (calc_history_reg(gshare_index,n) << (m-n))|calc_history_reg(index,m-n);
	this->predictions_values.prediction_event++;
    //printf("PRINTING \n");
	if((outcome == 't' && this->bk_h[final_index].counter < 2) | (outcome == 'n' && this->bk_h[final_index].counter > 1)) this->predictions_values.miss_prediction++;
	//printf("PRINTING \n");
    if(outcome == 't' && this->bk_h[final_index].counter != 3) this->bk_h[final_index].counter++;
	else if(outcome == 'n' && this->bk_h[final_index].counter != 0) this->bk_h[final_index].counter--;
	if(outcome == 't') this -> global_history_register = (this->global_history_register >> 1) | (1<<(n-1));
	else this->global_history_register = (this->global_history_register >> 1);
}

void hybrid_prediction(int i,predictor_table_t *this, unsigned long int addr, int M1, int M2, int N, int k,char outcome){
    this->predictions_values.prediction_event++;
    int index_ct      = calc_index(addr, k);
    int index_bimodal = calc_index(addr,M2);
    int index_gshare1  = calc_index(addr,M1);
    char b_outcome;
    char g_outcome;
    char h_outcome;

    // OBTAINING PREDICTION FROM BIMODAL 
    if(this->bk[index_bimodal].counter < 2) b_outcome = 'n'; 
    else if(this->bk[index_bimodal].counter > 1) b_outcome = 't';

    // OBTAINING PREDICTION FROM GSHARE
    unsigned long int n_part = bit_extraction(index_gshare1,M1,N);
	unsigned long int n_part_history_register = this->global_history_register;
	unsigned long int gshare_index = n_part_history_register ^ n_part;
	unsigned long int final_index = (calc_history_reg(gshare_index,N) << (M1-N))|calc_history_reg(index_gshare1,M1-N);

	if(this->bk_h[final_index].counter < 2) g_outcome = 'n';
    else if(this->bk_h[final_index].counter > 1) g_outcome = 't';
	
    // Updating Branch History Register
    if(outcome == 't') this -> global_history_register = (this->global_history_register >> 1) | (1<<(N-1));
	else this->global_history_register = (this->global_history_register >> 1);
    
    // Chooser index part

    if(this->bk_ct[index_ct].counter > 1){
        h_outcome = 'g';
        if(outcome == 't' && this->bk_h[final_index].counter != 3) this->bk_h[final_index].counter++;
	    else if(outcome == 'n' && this->bk_h[final_index].counter != 0) this->bk_h[final_index].counter--;
        if(g_outcome != outcome) this->predictions_values.miss_prediction++;

    } 
    else if(this->bk_ct[index_ct].counter < 2) {
        h_outcome = 'b';
        if(outcome == 't' && this->bk[index_bimodal].counter != 3) this->bk[index_bimodal].counter++;
	    else if(outcome == 'n' && this->bk[index_bimodal].counter != 0) this->bk[index_bimodal].counter--;
        if(b_outcome != outcome) this->predictions_values.miss_prediction++;
    
    }

    // Updating Chooser Table

    if(outcome == b_outcome && outcome != g_outcome && this->bk_ct[index_ct].counter != 0) this->bk_ct[index_ct].counter--;
    else if(outcome != b_outcome && outcome == g_outcome && this->bk_ct[index_ct].counter != 3) this->bk_ct[index_ct].counter++;
    //if(i < 70) {
      //  printf(" %d \t",i);
    //printf(" %x \t",addr);
       // printf("INDEX CT: %d \t",index_ct);
       // printf("CONTER VAL: %d \n",this->bk_ct[index_ct].counter);}
    




}



int main (int argc, char* argv[])
{
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file
	unsigned long int index;
	predictor_table_t normal_prediction;
    
    if (!(argc == 4 || argc == 5 || argc == 7))
    {
        printf("Error: Wrong number of inputs:%d\n", argc-1);
        exit(EXIT_FAILURE);
    }
    params.bp_name  = argv[1];
	//printf("bimodel \n");
    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if(strcmp(params.bp_name, "bimodal") == 0)              // Bimodal
    {
		//printf("BIMODAL \n");
        if(argc != 4)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M2       = strtoul(argv[2], NULL, 10);
        trace_file      = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
	
    }
    else if(strcmp(params.bp_name, "gshare") == 0)          // Gshare
    {
        if(argc != 5)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M1       = strtoul(argv[2], NULL, 10);
        params.N        = strtoul(argv[3], NULL, 10);
        trace_file      = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);

    }
    else if(strcmp(params.bp_name, "hybrid") == 0)          // Hybrid
    {
        if(argc != 7)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.K        = strtoul(argv[2], NULL, 10);
        params.M1       = strtoul(argv[3], NULL, 10);
        params.N        = strtoul(argv[4], NULL, 10);
        params.M2       = strtoul(argv[5], NULL, 10);
        trace_file      = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);

    }
    else
    {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }
    
    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if(FP == NULL)
    {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }
    
	//index = calc_index(addr,params.M2);
	if(strcmp(params.bp_name, "gshare") == 0) initialise_predictor_table(&normal_prediction,0,0,params.M1);
	else if(strcmp(params.bp_name, "bimodal") == 0) initialise_predictor_table(&normal_prediction,params.M2,0,0);
    else initialise_predictor_table(&normal_prediction,params.M2,params.K,params.M1);
    char str[2];
    int i;
    while(fscanf(FP, "%lx %s", &addr, str) != EOF)
    {
       	if(strcmp(params.bp_name, "bimodal") == 1) index = calc_index(addr,params.M1);
		else index = calc_index(addr,params.M2);
       	outcome = str[0];
		if(strcmp(params.bp_name, "gshare") == 0) gshare_prediction(&normal_prediction,index,params.M1,params.N,outcome);
		else if(strcmp(params.bp_name, "bimodal") == 0) prediction(&normal_prediction,index,outcome);
        else hybrid_prediction(i,&normal_prediction, addr, params.M1, params.M2, params.N, params.K,outcome);
        i++;
    }
	float miss_rate;
	miss_rate = ((float)(normal_prediction.predictions_values.miss_prediction))/((float)(normal_prediction.predictions_values.prediction_event)) * 100;
	
	printf("OUTPUT\n");
	printf(" number of predictions:		%d\n", normal_prediction.predictions_values.prediction_event);
	printf(" number of mispredictions:	%d\n", normal_prediction.predictions_values.miss_prediction);
	printf(" misprediction rate:		%.2f%\n",miss_rate);
	if(strcmp(params.bp_name, "bimodal") == 0) printf("FINAL BIMODAL CONTENTS\n");
	else if(strcmp(params.bp_name, "gshare") == 0) printf("FINAL GSHARE CONTENTS\n");
    else printf("FINAL CHOOSER CONTENTS\n");

	int j;
	int length;
	if(strcmp(params.bp_name, "bimodal") == 0){ 
        length = params.M2;
        //for(j = 0; j<calc_size(length); j++){
		  //  printf("%d \t", j);
		  //  printf("%d \n", normal_prediction.bk[j].counter);
	    //}
    }
	//else if(strcmp(params.bp_name, "gshare") == 0){
      //  length = params.M1;
        //for(j = 0; j<calc_size(length); j++){
		  //  printf("%d \t", j);
		   // printf("%d \n", normal_prediction.bk_h[j].counter);
        //}
    //}
    else if(strcmp(params.bp_name, "hybrid") == 0){
        length = params.K;
        for(j = 0; j<calc_size(length); j++){
		    printf("%d \t", j);
		    printf("%d \n", normal_prediction.bk_ct[j].counter);
        }
        printf("FINAL GSHARE CONTENTS\n");
        for(j = 0; j<calc_size(params.M1); j++){
		printf("%d \t", j);
		printf("%d \n", normal_prediction.bk_h[j].counter);}

        printf("FINAL BIMODAL CONTENTS\n");
        for(j = 0; j<calc_size(params.M2); j++){
		printf("%d \t", j);
		printf("%d \n", normal_prediction.bk[j].counter);}



    }
    return 0;
}
