#include <stdio.h>
#include "Eval.h"
#include <flint/fmpz.h>
#include <unistd.h>
int main() {
    Para *para = malloc(sizeof(*para));
    Data *data = malloc(sizeof(*data));
    SK *sk = malloc(sizeof(*sk));
    Gen_Data(para, data, sk);
    fmpz_t y;
    fmpz_init(y);
    int **index_S2;
    int **index_S3;
    int **index_C;
    int *index_sign_S2;
    int *index_sign_S3;
    int *index_sign_C;
    index_S2 = malloc(sizeof(int *) * 7);
    index_sign_S2 = malloc(sizeof(int) * 7);
    for (int i = 0; i < 7; i++) {
        index_S2[i] = malloc(sizeof(int) * 3);
    }
    index_S3 = malloc(sizeof(int *) * 12);
    index_sign_S3 = malloc(sizeof(int) * 12);
    for (int i = 0; i < 12; i++) {
        index_S3[i] = malloc(sizeof(int) * 3);
    }
    index_C = malloc(sizeof(int *) * 6);
    index_sign_C = malloc(sizeof(int) * 6);
    for (int i = 0; i < 6; i++) {
        index_C[i] = malloc(sizeof(int) * 3);
    }
    GetIndexTerm(index_S2, index_sign_S2, 2);
    GetIndexTerm(index_S3, index_sign_S3, 3);
    GetIndexTerm(index_C, index_sign_C, 4);
    clock_t S1_start, S1_finish, S2_start, S2_finish, S3_start, S3_finish, C_start, C_finish;
    double S1_time, S2_time, S3_time, C_time;
    //server 1
    ServerInput *input1 = malloc(sizeof(*input1));
    fmpz_mod_poly_init(input1->res, para->ctx);
    input1->para = para;
    input1->S = data->S1;
    input1->coeff = data->coeff;
    S1_start = clock();
    Eval1(input1);
    S1_finish = clock();
    S1_time = (double) (S1_finish - S1_start) / CLOCKS_PER_SEC;
    //server 2
    ServerInput *input2 = malloc(sizeof(*input2));
    fmpz_mod_poly_init(input2->res, para->ctx);
    input2->para = para;
    input2->S = data->S2;
    input2->coeff = data->coeff;
    input2->term_index = index_S2;
    input2->term_sign = index_sign_S2;
    S2_start = clock();
    Eval2(input2);
    S2_finish = clock();
    S2_time = (double) (S2_finish - S2_start) / CLOCKS_PER_SEC;
    //server 3
    ServerInput *input3 = malloc(sizeof(*input3));
    fmpz_mod_poly_init(input3->res, para->ctx);
    input3->para = para;
    input3->S = data->S3;
    input3->coeff = data->coeff;
    input3->term_index = index_S3;
    input3->term_sign = index_sign_S3;
    S3_start = clock();
    Eval3(input3);
    S3_finish = clock();
    S3_time = (double) (S3_finish - S3_start) / CLOCKS_PER_SEC;
    printf("all servers' running time: %f ms\n", (S1_time + S2_time + S3_time) * 1000);
    sleep(1);
    //client
    ClientInput *input = malloc(sizeof(*input));
    fmpz_init(input->y);
    fmpz_mod_poly_init(input->res1, para->ctx);
    fmpz_mod_poly_init(input->res2, para->ctx);
    fmpz_mod_poly_init(input->res3, para->ctx);
    input->para = para;
    input->sk = sk;
    input->data = data;
    input->term_index_S2 = index_S2;
    input->term_sign_S2 = index_sign_S2;
    input->term_index_S3 = index_S3;
    input->term_sign_S3 = index_sign_S3;
    input->term_index = index_C;
    input->term_sign = index_sign_C;
    C_start = clock();
    fmpz_mod_poly_set(input->res1, input1->res, para->ctx);
    fmpz_mod_poly_set(input->res2, input2->res, para->ctx);
    fmpz_mod_poly_set(input->res3, input3->res, para->ctx);
    Ver(input);
    C_finish = clock();
    C_time = (double) (C_finish - C_start) / CLOCKS_PER_SEC;
    printf("client running time : %f ms\n", C_time * 1000);
    return 0;
}
