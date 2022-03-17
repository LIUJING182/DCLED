#include <stdio.h>
#include "Eval.h"
#include <flint/fmpz.h>
#include <unistd.h>
#include "pthread.h"

int main() {
    Para *para = malloc(sizeof(*para));
    Data *data = malloc(sizeof(*data));
    SK *sk = malloc(sizeof(*sk));
    Gen_Data(para, data, sk);
    fmpz_t y;
    fmpz_init(y);
    fmpz_mod_poly_t res1, res2;
    fmpz_mod_poly_init(res1, para->ctx);
    fmpz_mod_poly_init(res2, para->ctx);
    clock_t S1_start, S1_finish, S2_start, S2_finish, C_start, C_finish;
    double S1_time, S2_time, C_time;
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
    sleep(1);
    //server 2
    ServerInput *input2 = malloc(sizeof(*input2));
    fmpz_mod_poly_init(input2->res, para->ctx);
    input2->para = para;
    input2->S = data->S2;
    input2->coeff = data->coeff;
    S2_start = clock();
    Eval2(input2);
    S2_finish = clock();
    S2_time = (double) (S2_finish - S2_start) / CLOCKS_PER_SEC;
    sleep(1);
    printf("all servers' running time: %f ms\n", (S1_time + S2_time) * 1000);
    //client
    ClientInput *input = malloc(sizeof(*input));
    fmpz_init(input->y);
    fmpz_mod_poly_init(input->res1, para->ctx);
    fmpz_mod_poly_init(input->res2, para->ctx);
    input->para = para;
    input->sk = sk;
    input->data = data;
    C_start = clock();
    fmpz_mod_poly_set(input->res1, input1->res, para->ctx);
    fmpz_mod_poly_set(input->res2, input2->res, para->ctx);
    Ver(input);
    C_finish = clock();
    C_time = (double) (C_finish - C_start) / CLOCKS_PER_SEC;
    printf("client's running time : %f ms\n", C_time * 1000);
    return 0;
}
