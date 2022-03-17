#ifndef DCLED_DEG2_EVAL_H
#define DCLED_DEG2_EVAL_H
#endif //DCLED_DEG2_EVAL_H

#include "tool.h"

typedef struct {
    fmpz_mod_poly_t res;
    Para *para;
    fmpz_mod_poly_t **S;
    fmpz_t *coeff;
} ServerInput;

typedef struct {
    fmpz_t y;
    fmpz_mod_poly_t res1;
    fmpz_mod_poly_t res2;
    Para *para;
    SK *sk;
    Data *data;
} ClientInput;

void *Eval1(void *input1) {
    ServerInput *input = (ServerInput *) input1;
    fmpz_mod_poly_t temp1;
    fmpz_mod_poly_t temp2;
    fmpz_mod_poly_init(temp1, input->para->ctx);
    fmpz_mod_poly_init(temp2, input->para->ctx);

    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            fmpz_mod_poly_mul(temp1, input->S[i][0], input->S[j][0], input->para->ctx);
            fmpz_mod_poly_mul(temp2, input->S[i][1], input->S[j][1], input->para->ctx);
            fmpz_mod_poly_sub(temp1, temp1, temp2, input->para->ctx);
            fmpz_mod_poly_scalar_mul_fmpz(temp1, temp1, input->coeff[i], input->para->ctx);
            fmpz_mod_poly_add(input->res, input->res, temp1, input->para->ctx);
        }
    }
    return NULL;
}

void *Eval2(void *input2) {
    ServerInput *input = (ServerInput *) input2;
    fmpz_mod_poly_t temp1;
    fmpz_mod_poly_t temp2;
    fmpz_mod_poly_init(temp1, input->para->ctx);
    fmpz_mod_poly_init(temp2, input->para->ctx);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            fmpz_mod_poly_mul(temp1, input->S[i][0], input->S[j][1], input->para->ctx);
            fmpz_mod_poly_mul(temp2, input->S[i][1], input->S[j][0], input->para->ctx);
            fmpz_mod_poly_add(temp1, temp1, temp2, input->para->ctx);
            fmpz_mod_poly_scalar_mul_fmpz(temp1, temp1, input->coeff[i], input->para->ctx);
            fmpz_mod_poly_add(input->res, input->res, temp1, input->para->ctx);
        }
    }
    return NULL;
}


void *Ver(void *input1) {
    ClientInput *input = (ClientInput *) input1;
    fmpz_t R1, R2, r1, r2, temp1, temp2, y1, y2;
    fmpz_init(R1);
    fmpz_init(R2);
    fmpz_init(r1);
    fmpz_init(r2);
    fmpz_init(temp1);
    fmpz_init(temp2);
    fmpz_init(y1);
    fmpz_init(y2);
    fmpz_mod_poly_evaluate_fmpz(r1, input->res1, input->sk->s1, input->para->ctx);
    fmpz_mod_poly_evaluate_fmpz(r2, input->res2, input->sk->s2, input->para->ctx);
    //R1
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            fmpz_mod_mul(temp1, input->data->r[i][0], input->data->r[j][0], input->para->ctx);
            fmpz_mod_mul(temp2, input->data->r[i][1], input->data->r[j][1], input->para->ctx);
            fmpz_mod_sub(temp1, temp1, temp2, input->para->ctx);
            fmpz_mod_mul(temp1, temp1, input->data->coeff[i], input->para->ctx);
            fmpz_mod_add(R1, R1, temp1, input->para->ctx);
        }
    }
    //R2
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            fmpz_mod_mul(temp1, input->data->r[i][2], input->data->r[j][3], input->para->ctx);
            fmpz_mod_mul(temp2, input->data->r[i][3], input->data->r[j][2], input->para->ctx);
            fmpz_mod_add(temp1, temp1, temp2, input->para->ctx);
            fmpz_mod_mul(temp1, temp1, input->data->coeff[i], input->para->ctx);
            fmpz_mod_add(R2, R2, temp1, input->para->ctx);
        }
    }
    //b
    fmpz_t b;
    fmpz_init(b);
    NativeEval(b, 2, input->para->num_data, input->data->coeff, input->data->b, input->para->ctx);
    if (fmpz_equal(r1, R1) == 1 && fmpz_equal(r2, R2) == 1) {
        fmpz_mod_poly_get_coeff_fmpz(y1, input->res1, 0, input->para->ctx);
        fmpz_mod_poly_get_coeff_fmpz(y2, input->res2, 0, input->para->ctx);
        fmpz_mod_add(input->y, y1, y2, input->para->ctx);
        fmpz_mod_add(input->y, input->y, b, input->para->ctx);
        printf("******************Verification Succeed!******************\n");
    } else {
        printf("******************Verification Failed!******************\n");
    }
}


