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
    for (int i = 0; i < m_rows; i ++) {    //seems like could replace with one forloop
        int k;   
        for (k = 0; k < m_cols/4*4; k += 4) {
            get_d[i*m_cols + k] = val;
            get_d[i*m_cols + k + 1] = val;
            get_d[i*m_cols + k + 2] = val;
            get_d[i*m_cols + k + 3] = val;
        }
        for (; k < m_cols; k++) {   //tail
            get_d[i*m_cols + k] = val;
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

    int rows_re = result->rows;
    int rows_a = mat1->rows;
    int rows_b = mat2->rows;

    if (cols_a != rows_b || rows_re != rows_a || cols_re != cols_b) {
        return -3;
    }
    double* data_a = mat1->data;
    double* data_b = mat2->data;
    double* data_re = result->data;

    fill_matrix(result, 0);

    #pragma omp parallel for
    for (int i = 0; i < rows_a; i++) {
        int k;
        #pragma omp parallel for 
        for (k = 0; k < cols_a/4*4; k += 4) {
            double one_row[cols_re];
            #pragma omp parallel for 
            for(int j = 0; j < cols_b; j++) {
                __m256d tmp_a = _mm256_loadu_pd(data_a + i*cols_a + k);   //split one row to each 4, 4, 4 ... items

                double b0 = data_b[k*cols_b + j];     //split b's col to each 4, 4, 4 ... items
                double b1 = data_b[(k+1)*cols_b + j];
                double b2 = data_b[(k+2)*cols_b + j];     //method 1
                double b3 = data_b[(k+3)*cols_b + j];
                __m256d tmp_b = _mm256_set_pd(b3, b2, b1, b0);
                __m256d tmp_sum = _mm256_mul_pd(tmp_a, tmp_b);
                one_row[j] = tmp_sum[0] + tmp_sum[1] + tmp_sum[2] + tmp_sum[3]; 
            }
            #pragma omp critical
            for (int j = 0; j < cols_b; j++) {
                data_re[i*cols_re + j] += one_row[j];
            }
        }
        #pragma omp parallel for 
        for (k = cols_a/4*4; k < cols_a; k++) {         //tail
            #pragma omp parallel for 
            for(int j = 0; j < cols_b; j++) {
                data_re[i*cols_re + j] += data_a[i*cols_a + k] * data_b[k*cols_b + j]; 
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
    /* A. FOR LOOP VERSION
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
    */

    // B. RECURSION VERSION 
    int cols_re = result->cols;
    int cols_a = mat->cols;

    int rows_re = result->rows;
    int rows_a = mat->rows;

    if (cols_re!= cols_a || rows_re != rows_a || pow < 0 || cols_a != rows_a) {
        return -3;
    }
    if (pow == 1) {         //if power = 1  =>  result.matrix = mat.matrix
        for (int i = 0; i < cols_re * rows_re; i++) {
            result->data[i] = mat->data[i];            //NEED SPEED UP
        }
        return 0;
    } else if (pow == 0) {  //if power = 0 -> result.data = mat.data
        for (int i = 0; i < rows_re; i++) {      //if power = 0  =>  identical matrix
            result->data[i*cols_re + i] = 1;     //identical matrix [[1,0],[0,1]];
        }
        return 0;
    } else {
        matrix *tmp_a = NULL;
        allocate_matrix(&tmp_a, rows_re, cols_re);
        
        matrix *tmp_b = NULL;
        allocate_matrix(&tmp_b, rows_re, cols_re);
        if ((pow & 1) == 0) {               //if pow is even  => result = (A^(n/2))^2
            pow_matrix(tmp_a, mat, pow/2);     
            mul_matrix(result, tmp_a, tmp_a);  
        } else if ((pow & 1) != 0) {        //if pow is odd   => result = A *(A^(n/2))^2
            pow_matrix(tmp_a, mat, pow/2);
            mul_matrix(tmp_b, tmp_a, tmp_a);
            mul_matrix(result, tmp_b, mat);    
        }
        deallocate_matrix(tmp_a);
        deallocate_matrix(tmp_b);
    }
    
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
