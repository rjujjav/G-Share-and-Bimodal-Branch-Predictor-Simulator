#include <stdlib.h>
#include <stdbool.h>
#ifndef SIM_BP_H
#define SIM_BP_H

typedef struct bp_params{
    unsigned long int K;
    unsigned long int M1;
    unsigned long int M2;
    unsigned long int N;
    char*             bp_name;
}bp_params;

typedef struct block{
	int counter;
} block_t;

typedef struct predict_stat{
	size_t miss_prediction;
	size_t prediction_event;
	size_t missprediction_rate; 
} predict_stat_t;

typedef struct predictor_table{
	int table_index;
	block_t *bk;
    block_t *bk_h;
    block_t *bk_ct;
	predict_stat_t predictions_values;
	unsigned long int global_history_register;
} predictor_table_t;
// Put additional data structures here as per your requirement

#endif
