#ifndef DCLED_DEG3_EVAL_H
#define DCLED_DEG3_EVAL_H
#endif //DCLED_DEG3_EVAL_H
#include "tool.h"
typedef struct {
    fmpz_mod_poly_t res;
    Para *para;
    fmpz_mod_poly_t **S;
    fmpz_t *coeff;
    int **term_index;
    int *term_sign;
} ServerInput;
typedef struct {
    fmpz_t y;
    fmpz_mod_poly_t res1;
    fmpz_mod_poly_t res2;
    fmpz_mod_poly_t res3;
    Para *para;
    SK *sk;
    Data *data;
    int **term_index_S2;
    int *term_sign_S2;
    int **term_index_S3;
    int *term_sign_S3;
    int **term_index;
    int *term_sign;
} ClientInput;

void *Eval1(void *input1) {
    ServerInput *input = (ServerInput *) input1;
    fmpz_mod_poly_t temp;
    fmpz_mod_poly_init(temp, input->para->ctx);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                fmpz_mod_poly_mul(temp, input->S[i][0], input->S[j][0], input->para->ctx);
                fmpz_mod_poly_mul(temp, temp, input->S[k][0], input->para->ctx);
                fmpz_mod_poly_scalar_mul_fmpz(temp, temp, input->coeff[i], input->para->ctx);
                fmpz_mod_poly_add(input->res, input->res, temp, input->para->ctx);
            }
        }
    }
    return NULL;
}

void *Eval2(void *input1) {
    ServerInput *input = (ServerInput *) input1;
    fmpz_mod_poly_t temp;
    fmpz_mod_poly_init(temp, input->para->ctx);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                for (int l = 0; l < 7; l++) {
                    fmpz_mod_poly_mul3_add(temp, temp, input->term_sign[l], input->S[i][input->term_index[l][0]],
                                           input->S[j][input->term_index[l][1]],
                                           input->S[k][input->term_index[l][2]], input->para->ctx);
                }
                fmpz_mod_poly_scalar_mul_fmpz(temp, temp, input->coeff[i], input->para->ctx);
                fmpz_mod_poly_add(input->res, input->res, temp, input->para->ctx);
                fmpz_mod_poly_init(temp, input->para->ctx);
            }
        }
    }
    return NULL;
}

void *Eval3(void *input1) {
    ServerInput *input = (ServerInput *) input1;
    fmpz_mod_poly_t temp;
    fmpz_mod_poly_init(temp, input->para->ctx);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                for (int l = 0; l < 12; l++) {
                    fmpz_mod_poly_mul3_add(temp, temp, input->term_sign[l],
                                           input->S[i][input->term_index[l][0]],
                                           input->S[j][input->term_index[l][1]],
                                           input->S[k][input->term_index[l][2]], input->para->ctx);
                }
                fmpz_mod_poly_scalar_mul_fmpz(temp, temp, input->coeff[i], input->para->ctx);
                fmpz_mod_poly_add(input->res, input->res, temp, input->para->ctx);
                fmpz_mod_poly_init(temp, input->para->ctx);
            }
        }
    }
    return NULL;
}

void *Ver(void *input1) {
    ClientInput *input = (ClientInput *) input1;
    fmpz_t R1, R2, R3, r1, r2, r3, temp1, temp2, y1, y2, y3;
    fmpz_init(R1);
    fmpz_init(R2);
    fmpz_init(R3);
    fmpz_init(r1);
    fmpz_init(r2);
    fmpz_init(r3);
    fmpz_init(temp1);
    fmpz_init(temp2);
    fmpz_init(y1);
    fmpz_init(y2);
    fmpz_init(y3);
    fmpz_mod_poly_evaluate_fmpz(r1, input->res1, input->sk->s1, input->para->ctx);
    fmpz_mod_poly_evaluate_fmpz(r2, input->res2, input->sk->s2, input->para->ctx);
    fmpz_mod_poly_evaluate_fmpz(r3, input->res3, input->sk->s3, input->para->ctx);
    //R1
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {


            for (int k = j; k < input->para->num_data; k++) {
                fmpz_mod_mul(temp1, input->data->r[i][0][0], input->data->r[j][0][0], input->para->ctx);
                fmpz_mod_mul(temp1, temp1, input->data->r[k][0][0], input->para->ctx);
                fmpz_mod_mul_fmpz(temp1, temp1, input->data->coeff[i], input->para->ctx);
                fmpz_mod_add(R1, R1, temp1, input->para->ctx);
            }
        }
    }
    //R2
    fmpz_set_ui(temp1, 0);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                for (int l = 0; l < 7; l++) {
                    fmpz_mod_mul3_add(temp1, temp1, input->term_sign_S2[l],
                                      input->data->r[i][1][input->term_index_S2[l][0]],
                                      input->data->r[j][1][input->term_index_S2[l][1]],
                                      input->data->r[k][1][input->term_index_S2[l][2]], input->para->ctx);
                }
                fmpz_mod_mul(temp1, temp1, input->data->coeff[i], input->para->ctx);
                fmpz_mod_add(R2, R2, temp1, input->para->ctx);
                fmpz_set_ui(temp1, 0);
            }
        }
    }
    //R3
    fmpz_set_ui(temp1, 0);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                for (int l = 0; l < 12; l++) {
                    fmpz_mod_mul3_add(temp1, temp1, input->term_sign_S3[l],
                                      input->data->r[i][2][input->term_index_S3[l][0]],
                                      input->data->r[j][2][input->term_index_S3[l][1]],
                                      input->data->r[k][2][input->term_index_S3[l][2]], input->para->ctx);
                }
                fmpz_mod_mul(temp1, temp1, input->data->coeff[i], input->para->ctx);
                fmpz_mod_add(R3, R3, temp1, input->para->ctx);
                fmpz_set_ui(temp1, 0);
            }
        }
    }
    //b
    fmpz_t b;
    fmpz_init(b);
    fmpz_set_ui(temp1, 0);
    for (int i = 0; i < input->para->num_data; i++) {
        for (int j = i; j < input->para->num_data; j++) {
            for (int k = j; k < input->para->num_data; k++) {
                for (int l = 0; l < 6; l++) {
                    fmpz_mod_mul3_add(temp1, temp1, input->term_sign[l],
                                      input->data->a[i][input->term_index[l][0]],
                                      input->data->a[j][input->term_index[l][1]],
                                      input->data->a[k][input->term_index[l][2]], input->para->ctx);
                }
                fmpz_mod_mul(temp1, temp1, input->data->coeff[i], input->para->ctx);
                fmpz_mod_add(b, b, temp1, input->para->ctx);
                fmpz_set_ui(temp1, 0);
            }
        }
    }
    if (fmpz_equal(r1, R1) == 1 && fmpz_equal(r2, R2) == 1 && fmpz_equal(r3, R3) == 1) {
        fmpz_mod_poly_get_coeff_fmpz(y1, input->res1, 0, input->para->ctx);
        fmpz_mod_poly_get_coeff_fmpz(y2, input->res2, 0, input->para->ctx);
        fmpz_mod_poly_get_coeff_fmpz(y3, input->res3, 0, input->para->ctx);
        fmpz_mod_add(input->y, y1, y2, input->para->ctx);
        fmpz_mod_add(input->y, input->y, y3, input->para->ctx);
        fmpz_mod_add(input->y, input->y, b, input->para->ctx);
        printf("******************Verification Succeed!******************\n");
    } else {
        printf("******************Verification Failed!******************\n");
    }
}

