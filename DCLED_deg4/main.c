#include <stdio.h>
#include "Eval.h"
#include <flint/fmpz.h>

int main() {
    Para *para = malloc(sizeof(*para));
    Data *data = malloc(sizeof(*data));
    SK *sk = malloc(sizeof(*sk));
    Gen_Data(para, data, sk);
    fmpz_t y;
    fmpz_init(y);
    int **index_S2;
    int **index_S3;
    int **index_S4;
    int **index_C;
    int *index_sign_S2;
    int *index_sign_S3;
    int *index_sign_S4;
    int *index_sign_C;
    index_S2 = malloc(sizeof(int *) * 45);
    index_sign_S2 = malloc(sizeof(int) * 45);
    for (int i = 0; i < 45; i++) {
        index_S2[i] = malloc(sizeof(int) * 4);
    }
    index_S3 = malloc(sizeof(int *) * 300);
    index_sign_S3 = malloc(sizeof(int) * 300);
    for (int i = 0; i < 300; i++) {
        index_S3[i] = malloc(sizeof(int) * 4);
    }
    index_S4 = malloc(sizeof(int *) * 360);
    index_sign_S4 = malloc(sizeof(int) * 360);
    for (int i = 0; i < 360; i++) {
        index_S4[i] = malloc(sizeof(int) * 4);
    }
    index_C = malloc(sizeof(int *) * 23);
    index_sign_C = malloc(sizeof(int) * 23);
    for (int i = 0; i < 23; i++) {
        index_C[i] = malloc(sizeof(int) * 4);
    }
    GetIndexTerm(index_S2, index_sign_S2, 2);
    GetIndexTerm(index_S3, index_sign_S3, 3);
    GetIndexTerm(index_S4, index_sign_S4, 4);
    GetIndexTerm(index_C, index_sign_C, 5);
    clock_t S_start, S_finish, C_start, C_finish;
    double S_time, C_time;
    //server 1
    ServerInput *input1 = malloc(sizeof(*input1));
    fmpz_mod_poly_init(input1->res, para->ctx);
    input1->para = para;
    input1->S = data->S1;
    input1->coeff = data->coeff;
    //server 2
    ServerInput *input2 = malloc(sizeof(*input2));
    fmpz_mod_poly_init(input2->res, para->ctx);
    input2->para = para;
    input2->S = data->S2;
    input2->coeff = data->coeff;
    input2->term_index = index_S2;
    input2->term_sign = index_sign_S2;
    //server 3
    ServerInput *input3 = malloc(sizeof(*input3));
    fmpz_mod_poly_init(input3->res, para->ctx);
    input3->para = para;
    input3->S = data->S3;
    input3->coeff = data->coeff;
    input3->term_index = index_S3;
    input3->term_sign = index_sign_S3;
    //server 4
    ServerInput *input4 = malloc(sizeof(*input4));
    fmpz_mod_poly_init(input4->res, para->ctx);
    input4->para = para;
    input4->S = data->S4;
    input4->coeff = data->coeff;
    input4->term_index = index_S4;
    input4->term_sign = index_sign_S4;
    S_start = clock();
    Eval1(input1);
    Eval2(input2);
    Eval3(input3);
    Eval4(input4);
    S_finish = clock();
    S_time = (double) (S_finish - S_start) / CLOCKS_PER_SEC;
    printf("all servers' running time: %f ms\n", S_time * 1000);
    sleep(1);
    //client
    ClientInput *input = malloc(sizeof(*input));
    fmpz_init(input->y);
    fmpz_mod_poly_init(input->res1, para->ctx);
    fmpz_mod_poly_init(input->res2, para->ctx);
    fmpz_mod_poly_init(input->res3, para->ctx);
    fmpz_mod_poly_init(input->res4, para->ctx);
    input->para = para;
    input->sk = sk;
    input->data = data;
    input->term_index_S2 = index_S2;
    input->term_sign_S2 = index_sign_S2;
    input->term_index_S3 = index_S3;
    input->term_sign_S3 = index_sign_S3;
    input->term_index_S4 = index_S4;
    input->term_sign_S4 = index_sign_S4;
    input->term_index = index_C;
    input->term_sign = index_sign_C;
    C_start = clock();
    fmpz_mod_poly_set(input->res1, input1->res, para->ctx);
    fmpz_mod_poly_set(input->res2, input2->res, para->ctx);
    fmpz_mod_poly_set(input->res3, input3->res, para->ctx);
    fmpz_mod_poly_set(input->res4, input4->res, para->ctx);
    Ver(input);
    C_finish = clock();
    C_time = (double) (C_finish - C_start) / CLOCKS_PER_SEC;
    printf("client's running time : %f ms\n", C_time * 1000);
    sleep(1);
    return 0;
}
