#ifndef DCLED_DEG3_TOOL_H
#define DCLED_DEG3_TOOL_H
#endif //DCLED_DEG3_TOOL_H

#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>

typedef struct {
    fmpz_t *x;
    fmpz_t **a;
    fmpz_t ***r;
    fmpz_mod_poly_t **S1;
    fmpz_mod_poly_t **S2;
    fmpz_mod_poly_t **S3;
    fmpz_t *coeff;
} Data;
typedef struct {
    fmpz_t s1;
    fmpz_t s2;
    fmpz_t s3;
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
    //init
    fmpz_mod_ctx_init(para->ctx, para->q);
    flint_rand_t state;
    flint_randinit(state);
    data->x = malloc(sizeof(fmpz_t) * para->num_data);
    data->coeff = malloc(sizeof(fmpz_t) * para->num_data);
    data->a = malloc(sizeof(fmpz_t *) * para->num_data);
    data->r = malloc(sizeof(fmpz_t **) * para->num_data);
    data->S1 = malloc(sizeof(fmpz_mod_poly_t *) * para->num_data);
    data->S2 = malloc(sizeof(fmpz_mod_poly_t *) * para->num_data);
    data->S3 = malloc(sizeof(fmpz_mod_poly_t *) * para->num_data);
    for (int i = 0; i < para->num_data; i++) {
        fmpz_init(data->x[i]);
        fmpz_init(data->coeff[i]);
        fmpz_randm(data->x[i], state, para->q);
        fmpz_randm(data->coeff[i], state, para->q);
        data->S1[i] = malloc(sizeof(fmpz_mod_poly_t) * 3);
        data->S2[i] = malloc(sizeof(fmpz_mod_poly_t) * 3);
        data->S3[i] = malloc(sizeof(fmpz_mod_poly_t) * 3);
        for (int j = 0; j < 3; j++) {
            fmpz_mod_poly_init(data->S1[i][j], para->ctx);
            fmpz_mod_poly_init(data->S2[i][j], para->ctx);
            fmpz_mod_poly_init(data->S3[i][j], para->ctx);
        }
    }
    for (int i = 0; i < para->num_data; i++) {
        data->a[i] = malloc(sizeof(fmpz_t) * 3);
        for (int j = 0; j < 3; j++) {
            fmpz_init(data->a[i][j]);
            fmpz_randm(data->a[i][j], state, para->q);
        }
        data->r[i] = malloc(sizeof(fmpz_t *) * 3);
        for (int j = 0; j < 3; j++) {
            data->r[i][j] = malloc(sizeof(fmpz_t) * 3);
            for (int k = 0; k < 3; k++) {
                fmpz_init(data->r[i][j][k]);
                fmpz_randm(data->r[i][j][k], state, para->q);
            }
        }
    }
    fmpz_init(sk->s1);
    fmpz_init(sk->s2);
    fmpz_init(sk->s3);
    fmpz_randm(sk->s1, state, para->q);
    fmpz_randm(sk->s2, state, para->q);
    fmpz_randm(sk->s3, state, para->q);
    //gen cipher of server S1 S2 S3
    fmpz_t inv_s1, inv_s2, inv_s3;
    fmpz_init(inv_s1);
    fmpz_init(inv_s2);
    fmpz_init(inv_s3);
    fmpz_invmod(inv_s1, sk->s1, para->q);
    fmpz_invmod(inv_s2, sk->s2, para->q);
    fmpz_invmod(inv_s3, sk->s3, para->q);
    fmpz_t temp;
    fmpz_init(temp);
    for (int i = 0; i < para->num_data; i++) {
        //s1
        //xi-ai1
        fmpz_mod_sub(temp, data->x[i], data->a[i][0], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][0], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][0][0], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s1, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][0], 1, temp, para->ctx);
        //ai2
        fmpz_set(temp, data->a[i][1]);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][1], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][0][1], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s1, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][1], 1, temp, para->ctx);
        //ai3
        fmpz_set(temp, data->a[i][2]);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][2], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][0][2], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s1, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S1[i][2], 1, temp, para->ctx);
        //s2
        //ai1
        fmpz_set(temp, data->a[i][0]);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][0], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][1][0], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s2, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][0], 1, temp, para->ctx);
        //xi-ai2
        fmpz_mod_sub(temp, data->x[i], data->a[i][1], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][1], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][1][1], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s2, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][1], 1, temp, para->ctx);
        //ai3
        fmpz_set(temp, data->a[i][2]);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][2], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][1][2], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s2, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S2[i][2], 1, temp, para->ctx);
        //s3
        //ai1
        fmpz_set(temp, data->a[i][0]);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][0], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][2][0], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s3, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][0], 1, temp, para->ctx);
        //ai2
        fmpz_set(temp, data->a[i][1]);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][1], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][2][1], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s3, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][1], 1, temp, para->ctx);
        //xi-ai3
        fmpz_mod_sub(temp, data->x[i], data->a[i][0], para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][2], 0, temp, para->ctx);
        fmpz_mod_sub(temp, data->r[i][2][2], temp, para->ctx);
        fmpz_mod_mul(temp, temp, inv_s3, para->ctx);
        fmpz_mod_poly_set_coeff_fmpz(data->S3[i][2], 1, temp, para->ctx);
    }
}

void GetIndexTerm(int **index_term, int *index_sign, int index) {
    FILE *fp;
    fp = fopen("3Severs.txt", "r");
    int a[78];
    for (int i = 0; i < 78; i++) {
        fscanf(fp, "%d", &a[i]);
    }
    fclose(fp);
    fp = fopen("3Servers-sign.txt", "r");
    int b[25];
    for (int i = 0; i < 25; i++) {
        fscanf(fp, "%d", &b[i]);
    }
    fclose(fp);
    int index_a;
    if (index == 2) {
        index_a = 3;
        for (int i = 0; i < 7; i++) {
            index_sign[i] = b[i];
            for (int j = 0; j < 3; j++) {
                index_term[i][j] = a[index_a];
                index_a = index_a + 1;
            }
        }
    }
    if (index == 3) {
        index_a = 24;
        for (int i = 0; i < 12; i++) {
            index_sign[i] = b[i + 7];
            for (int j = 0; j < 3; j++) {
                index_term[i][j] = a[index_a];
                index_a = index_a + 1;
            }
        }
    }
    if (index == 4) {
        index_a = 60;
        for (int i = 0; i < 6; i++) {
            index_sign[i] = b[i + 19];
            for (int j = 0; j < 3; j++) {
                index_term[i][j] = a[index_a];
                index_a = index_a + 1;
            }
        }
    }
}

void fmpz_mod_poly_mul3_add(fmpz_mod_poly_t res, fmpz_mod_poly_t add1, int sign, fmpz_mod_poly_t poly1,
                            fmpz_mod_poly_t poly2, fmpz_mod_poly_t poly3, fmpz_mod_ctx_t ctx) {
    fmpz_mod_poly_t temp;
    fmpz_mod_poly_init(temp, ctx);
    fmpz_mod_poly_mul(temp, poly1, poly2, ctx);
    fmpz_mod_poly_mul(temp, temp, poly3, ctx);
    if (sign == 1) {
        fmpz_mod_poly_add(res, add1, temp, ctx);
    } else {
        fmpz_mod_poly_sub(res, add1, temp, ctx);
    }
    fmpz_mod_poly_clear(temp, ctx);
}

void fmpz_mod_mul3_add(fmpz_t res, fmpz_t add1, int sign, fmpz_t poly1,
                       fmpz_t poly2, fmpz_t poly3, fmpz_mod_ctx_t ctx) {
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_mod_mul(temp, poly1, poly2, ctx);
    fmpz_mod_mul(temp, temp, poly3, ctx);
    if (sign == 1) {
        fmpz_mod_add(res, add1, temp, ctx);
    } else {
        fmpz_mod_sub(res, add1, temp, ctx);
    }
    fmpz_clear(temp);
}