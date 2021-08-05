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
    if (!new_mat->data) {
        return -2;
    }
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
    if (cols <= 0 || rows <= 0 ) {      //MAY NEED REVISE LATER?????
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
    //parent == null && has no slice (ref_cnt == 1)  ->  free `mat->data`
    if ((!mat->parent) && (mat->ref_cnt == 1)) {
        free(mat->data);
        free(mat);
    } else if ((!mat->parent) && (mat->ref_cnt > 1)) {
        mat->ref_cnt--;
    }
    if (!mat->parent) {
        return;
    }
    //`mat` is only child of its parent  &  parent's parent == null      -> free `mat->parent->data`
    else if (mat->parent && mat->parent->ref_cnt == 1) {
        free(mat->parent->data);
        free(mat);
        free(mat->parent);
    }
    else if (mat->parent && mat->parent->ref_cnt > 1) {
        mat->parent->ref_cnt--;
        free(mat);
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
    return get_d[col + row*m_cols];
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

    #pragma omp parallel for
    for (int i = 0; i < m_rows * m_cols; i ++) {    //seems like could replace with one forloop
        get_d[i] = val;
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

    #pragma omp parallel for
    for (int i = 0; i < rows_a; i ++) {
        int k;
        for (k = 0; k < cols_a/4*4; k += 4) {
            __m256d tmp_a = _mm256_loadu_pd(data_a + i*cols_a + k);
            __m256d tmp_b = _mm256_loadu_pd(data_b + i*cols_b + k);
            __m256d sums = _mm256_add_pd(tmp_a, tmp_b);

            
            data_re[i*cols_a + k] = sums[0];        //haven't decide if use storeu
            data_re[i*cols_a + k + 1] = sums[1];
            data_re[i*cols_a + k + 2] = sums[2];
            data_re[i*cols_a + k + 3] = sums[3];
        }
        for (; k < cols_a; k++) {         //tail
            data_re[i*cols_a + k] = data_a[i*cols_a + k] + data_b[i*cols_a + k];
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
    int rows_a = mat1->rows;
    int rows_b = mat2->rows;
    if (cols_a != mat2->rows || result->rows != rows_a || cols_re != cols_b) {
        return -3;
    }

    double* data_a = mat1->data;
    double* data_b = mat2->data;
    double* data_re = result->data;

    //fill_matrix(result, 0);
    matrix *trans;                            //transpose cite: https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
    allocate_matrix(&trans, cols_b, cols_a);  //error check?
    double* data_tran = trans->data;
    #pragma omp parallel for      
    for (int n = 0; n < cols_b * rows_b; n++) {
        data_tran[n] = data_b[cols_b*(n%cols_a) + n/cols_a];
    }

    #pragma omp parallel for if (rows_a >= 256 && cols_a >= 4)
    for (int i = 0; i < rows_a; i++) {
        #pragma omp parallel for if (cols_b >= 256 && cols_a >= 4)
        for(int j = 0; j < cols_b; j++) {
            int k;
            double dot_sum = 0.0;
            __m256d tmp_a;
            __m256d tmp_b;
            __m256d tmp_sum = _mm256_set1_pd(0);
            for (k = 0; k < cols_a/4*4; k += 4) {
                tmp_a = _mm256_loadu_pd(data_a + i*cols_a + k);   //split one row to 4, 4, 4 ... each
                tmp_b = _mm256_loadu_pd(data_tran + j*cols_a + k); 
                tmp_sum = _mm256_fmadd_pd(tmp_a, tmp_b, tmp_sum);
            }
            dot_sum += tmp_sum[0] + tmp_sum[1] + tmp_sum[2] + tmp_sum[3];
            for (k = cols_a/4*4; k < cols_a; k ++) {
                dot_sum += data_a[i*cols_a + k] * data_tran[j*cols_a + k];
            }
            data_re[i*cols_re + j] = dot_sum;
        }
    }
    deallocate_matrix(trans);
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    /**Citation: https://xlinux.nist.gov/dads/HTML/repeatedSquaring.html */
    int cols_re = result->cols;
    int cols_a = mat->cols;

    int rows_re = result->rows;
    int rows_a = mat->rows;
    if (cols_re != cols_a || rows_re != rows_a || pow < 0 || cols_a != rows_a) {
        return -3;
    }
    int pow_bianry = pow_to_binary(pow);
    matrix  *cache = NULL;
    matrix *tmp = NULL;
    matrix *cur = NULL;
    int a = allocate_matrix(&tmp, rows_re, cols_re);
    int b = allocate_matrix(&cache, rows_re, cols_re);
    int c = allocate_matrix(&cur, rows_re, cols_re);
    if (a || b || c) {
        return -2;
    }
    fill_matrix(result, 0);
    for (int i = 0; i < rows_re; i++) {
        result->data[i*cols_re + i] = 1;            // tmp =  identical matrix [[1,0],[0,1]];
    }
    
    int position = 0;

    while(pow_bianry) {
        if (!position) {
            mul_matrix(tmp, result, mat);
        } else {
            mul_matrix(tmp, cache, cache);
        }

        if (pow_bianry % 10 && !position) {
            mul_matrix(cur, tmp, result);
        } else if (pow_bianry % 10 && position) {
            mul_matrix(cur, result, tmp);
        }
  
        for (int j = 0; j < rows_re * cols_re; j++) {   
            if (pow_bianry % 10) {
                result->data[j] = cur->data[j];
            }
            cache->data[j] = tmp->data[j];
            
        }   
        position++;
        pow_bianry /= 10;
    }

    deallocate_matrix(cache);
    deallocate_matrix(tmp);
    deallocate_matrix(cur);

    return 0;
}

int pow_to_binary(int pow) {
    int binary = 0;
    int position = 1;
    while (pow) {
        if (pow % 2) {
            binary += 1 * position;
        } 
        position *= 10;
        pow /= 2;
    }
    return binary;
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

    #pragma omp parallel for
    for (int i = 0; i < m_rows; i ++) {
        int k;   
        for (k = 0; k < m_cols/4*4; k += 4) {
            __m256d neg_one = _mm256_set1_pd(-1);
            __m256d source = _mm256_loadu_pd(get_d + i*m_rows + k);
            __m256d neg_vec = _mm256_mul_pd(source, neg_one);
            _mm256_storeu_pd(re_d + i*m_cols + k, neg_vec);
        }
        for (; k < m_cols; k++) {   //tail
            re_d[i*m_cols + k] = get_d[i*m_rows + k]*(-1);
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

    double *re_d = result->data;

    #pragma omp parallel for
    for (int i = 0; i < m_rows; i++) {
        #pragma omp parallel for
        for (int j = 0; j < m_cols; j++) {
            re_d[i*m_cols + j] = fabs(get_d[i*m_cols + j]);
        }
    }
    return 0;
}
