#ifndef DCLED_DEG2_TOOL_H
#define DCLED_DEG2_TOOL_H
#endif //DCLED_DEG2_TOOL_H

#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

typedef struct {
    fmpz_t *x;
    fmpz_t *a;
    fmpz_t *b;
    fmpz_t **r;
    fmpz_mod_poly_t **S1;
    fmpz_mod_poly_t **S2;
    fmpz_t *coeff;
    int **arr;
} Data;
typedef struct {
    fmpz_t s1;
    fmpz_t s2;
} SK;
typedef struct {
    int num_data;
    fmpz_t q;
    fmpz_mod_ctx_t ctx;
} Para;

void Gen_Data(Para *para, Data *data, SK *sk) {
    //set num of data & msg space Zq
    printf("Enter n:\n");
    scanf("%d", &para->num_data);
    fmpz_init(para->q);
    fmpz_set_str(para->q, "340282366920938463463374607431768211507", 10);
    fmpz_mod_ctx_init(para->ctx, para->q);
    flint_rand_t state;
    flint_randinit(state);
    data->x = malloc(sizeof(fmpz_t) * para->num_data);
    data->coeff = malloc(sizeof(fmpz_t) * para->num_data);
    data->a = malloc(sizeof(fmpz_t) * para->num_data);
    data->b = malloc(sizeof(fmpz_t) * para->num_data);
    data->r = malloc(sizeof(fmpz_t *) * para->num_data);
    data->S1 = malloc(sizeof(fmpz_mod_poly_t *) * para->num_data);
    data->S2 = malloc(sizeof(fmpz_mod_poly_t *) * para->num_data);
    for (int i = 0; i < para->num_data; i++) {
        fmpz_init(data->x[i]);
        fmpz_init(data->coeff[i]);
        fmpz_init(data->a[i]);
        fmpz_init(data->b[i]);
        data->r[i] = malloc(sizeof(fmpz_t) * 4);
        fmpz_randm(data->x[i], state, para->q);
        fmpz_randm(data->coeff[i], state, para->q);
        fmpz_randm(data->a[i], state, para->q);
        fmpz_randm(data->b[i], state, para->q);
        for (int j = 0; j < 4; j++) {
            fmpz_init(data->r[i][j]);
            fmpz_randm(data->r[i][j], state, para->q);
        }
        data->S1[i] = malloc(sizeof(fmpz_mod_poly_t) * 2);
        data->S2[i] = malloc(sizeof(fmpz_mod_poly_t) * 2);
        for (int j = 0; j < 2; j++) {
            fmpz_mod_poly_init(data->S1[i][j], para->ctx);
            fmpz_mod_poly_init(data->S2[i][j], para->ctx);
        }
    }
    fmpz_init(sk->s1);
    fmpz_init(sk->s2);
    fmpz_randm(sk->s1, state, para->q);
    fmpz_randm(sk->s2, state, para->q);
    //gen cipher of server S1 S2
    fmpz_t inv_s1, inv_s2;
    fmpz_init(inv_s1);
    fmpz_init(inv_s2);
    fmpz_invmod(inv_s1, sk->s1, para->q);
    fmpz_invmod(inv_s2, sk->s2, para->q);
    fmpz_t temp;
    fmpz_init(temp);
    for (int i = 0; i < para->num_data; i++) {
        //x-a
        fmpz_mod_sub(temp, data->x[i], data->a[i], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][0], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][0], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s1, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][0], 1, temp, para->ctx);
        //a-b
        fmpz_mod_sub(temp, data->a[i], data->b[i], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][1], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][1], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s1, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][1], 1, temp, para->ctx);
        //m-b
        fmpz_mod_sub(temp, data->x[i], data->b[i], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][0], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][2], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s2, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][0], 1, temp, para->ctx);
        //a
        fmpz_set(temp, data->a[i]);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][1], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][3], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s2, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][1], 1, temp, para->ctx);
    }
}

void NativeEval_f(fmpz_t y, int d, int num_data, int loop, int beg_ind, int *ind_var, fmpz_t *coeff, fmpz_t *X,
                  fmpz_mod_ctx_t ctx) {
    if (loop == d) {
        fmpz_t temp;
        fmpz_init_set_ui(temp, 1);
        for (int i = 0; i < d; i++) {
            fmpz_mod_mul(temp, temp, X[ind_var[i]], ctx);
        }
        fmpz_mod_mul(temp, temp, coeff[ind_var[0]], ctx);
        fmpz_mod_add(y, y, temp, ctx);

    } else {
        loop = loop + 1;
        for (int i = beg_ind; i < num_data; i++) {
            ind_var[loop - 1] = i;
            NativeEval_f(y, d, num_data, loop, i, ind_var, coeff, X, ctx);
        }
    }
}

void NativeEval(fmpz_t y, int d, int num_data, fmpz_t *coeff, fmpz_t *X, fmpz_mod_ctx_t ctx) {
    fmpz_t temp;
    fmpz_init(temp);
    for (int i = 1; i < d + 1; i++) {
        int *ind_var = (int *) malloc(sizeof(int) * i);
        NativeEval_f(temp, i, num_data, 0, 0, ind_var, coeff, X, ctx);
        fmpz_mod_add(y, y, temp, ctx);
        fmpz_set_ui(temp, 0);
    }
}