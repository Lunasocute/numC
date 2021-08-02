#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails. Remember to set the error messages in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows <= 0 || cols <= 0) {
        return -1;
    }
    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (!new_mat) {
        return -2;
    }
    new_mat->rows = rows;
    new_mat->cols = cols;
    new_mat->ref_cnt = 1;
    new_mat->parent = NULL;
    
    new_mat->data = (double *) calloc(rows*cols, sizeof(double));
    *mat = new_mat;
    return 0;
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Remember to set the error messages in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (cols <= 0 || rows <= 0 ) {    //MAY NEED REVISE LATER
        return -1;
    }
    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (!new_mat) {
        return -2;
    }
    new_mat->rows = rows;
    new_mat->cols = cols;
    new_mat->ref_cnt = 1;
    new_mat->parent = from;
    new_mat->data = from->data + offset;  //pointer + offset
    *mat = new_mat;

    (from->ref_cnt)++;
    return 0;
}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or that you free `mat->parent->data` if `mat` is the last existing slice of its parent matrix and its parent matrix has no other references
 * (including itself). You cannot assume that mat is not NULL.
 */

/** check ref_cnt */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (!mat) {
        return;
    }
    // parent == null && has no slice (ref_cnt == 1)  ->  free `mat->data`
    if ((!mat->parent) && mat->ref_cnt == 1) {
        free(mat->data);
    }
    //`mat` is only child of its parent  &  parent's parent == null      -> free `mat->parent->data`
    if (!mat->parent) {
        return;
    }
    else if ((!mat->parent->parent) && mat->parent->ref_cnt == 2) {
        free(mat->parent->data);
    }
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    double* get_d = mat->data;
    int m_cols = mat->cols;
    
    double get_val = get_d[col + row*m_cols];
    return get_val;
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */

    int m_cols = mat->cols;
    double* get_d = mat->data;
    
    get_d[col + row*m_cols] = val;
}

/*
 * Sets all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    double* get_d = mat->data;
    int m_rows = mat->rows;
    int m_cols = mat->cols;

    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            get_d[i*m_cols + j] = val;
        }
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int cols_re = result->cols;
    int cols_a = mat1->cols;
    int cols_b = mat2->cols;

    int rows_re = result->rows;
    int rows_a = mat1->rows;
    int rows_b = mat2->rows;

    if (cols_a != cols_b || cols_re!= cols_b || rows_a != rows_b || rows_re != rows_b) {
        return -3;
    }

    double* data_a = mat1->data;
    double* data_b = mat2->data;
    double* data_re = result->data;


    for (int i = 0; i < rows_a; i++) {
        for (int j = 0; j < cols_a; j++) {
            data_re[i*cols_a + j] = data_a[i*cols_a + j] + data_b[i*cols_a + j];
        }
    }
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int cols_re = result->cols;
    int cols_a = mat1->cols;
    int cols_b = mat2->cols;
    int rows_re = result->rows;
    int rows_a = mat1->rows;
    int rows_b = mat2->rows;

    if (cols_a != rows_b || rows_re != rows_a || cols_re != cols_b) {
        return -3;
    }
    double* data_a = mat1->data;
    double* data_b = mat2->data;
    double* data_re = result->data;

    for (int i = 0; i < rows_re; i++) {
        for (int j = 0; j < cols_re; j++) {
            data_re[i*cols_re + j] = 0;
            for (int k = 0; k < rows_b; k++) {
                data_re[i*cols_re + j] += data_a[i*cols_a+ k] * data_b[j + k*cols_b];
            }
        }
    }

    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    int cols_re = result->cols;
    int cols_a = mat->cols;

    int rows_re = result->rows;
    int rows_a = mat->rows;
    if (cols_re != cols_a || rows_re != rows_a || pow < 0 || cols_a != rows_a) {
        return -3;
    }

    matrix *tmp = NULL;
    allocate_matrix(&tmp, rows_re, cols_re);
    fill_matrix(tmp, 0);
    for (int i = 0; i < rows_re; i++) {
        tmp->data[i*cols_re + i] = 1;            // tmp =  identical matrix [[1,0],[0,1]];
    }

    for (int i = 0; i < pow; i++) {
        mul_matrix(result, tmp, mat);
        for (int j = 0; j < rows_re * cols_re; j++) {    //NEED IMPROVE SPEED LATER
            tmp->data[j] = result->data[j];
        }   
    }
    deallocate_matrix(tmp);
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    double* get_d = mat->data;
    int m_rows = mat->rows;
    int m_cols = mat->cols;

    double* re_d = result->data;
    int re_rows = result->rows;
    int re_cols = result->cols;
    if (re_rows != m_rows || re_cols != m_cols){
        return -3;
    }

    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            re_d[i * m_rows + j] = get_d[i * m_cols + j] * (-1);
        }
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    double* get_d = mat->data;
    int m_rows = mat->rows;
    int m_cols = mat->cols;

    matrix *neg_ma = NULL;
    allocate_matrix(&neg_ma, m_rows, m_cols);
    double *re_d = result->data;


    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            re_d[i*m_cols + j] = abs(get_d[i*m_cols + j]);
        }
    }
    return 0;
}
