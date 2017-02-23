#include "DS.h"

static stack_datastruct stackinstant = { NULL, NULL, 0 };

void init_stack
    (
    void
    )
{
    stackinstant.stack_starting_address = _aligned_malloc(STACK_TOTAL_SIZE, 0x20);
    if (stackinstant.stack_starting_address == NULL)
        {
        printf("failed creating memory stack\n");
        exit(EXIT_FAILURE);
        }
    stackinstant.stack_current_address = (int8_t*)stackinstant.stack_starting_address;
    stackinstant.stack_current_alloc_size = 0;
}

void free_stack
    (
    void
    )
{
    _aligned_free(stackinstant.stack_starting_address);
}

void reset_stack_ptr_to_initial_position
    (
    void
    )
{
    stackinstant.stack_current_address = (char*)stackinstant.stack_starting_address;
    stackinstant.stack_current_alloc_size = 0;
}

void* alloc_from_stack
    (
    uint32_t len
    )
{
    void* ptr = NULL;
    if (len <= 0)
        {
        len = 0x20;
        }
    uint32_t aligned_len = (len + 0xF) & (~0xF);
    stackinstant.stack_current_alloc_size += aligned_len;
    if (stackinstant.stack_current_alloc_size >= STACK_TOTAL_SIZE)
        {
        printf("failed allocating memory from stack anymore\n");
        _aligned_free(stackinstant.stack_starting_address);
        exit(EXIT_FAILURE);
        }
    ptr = stackinstant.stack_current_address;
    stackinstant.stack_current_address += aligned_len;
    /* C99: all zero bits means 0 for fixed points, 0.0 for floating points */
    memset(ptr, 0, len);
    return ptr;
}

void partial_free_from_stack
    (
    uint32_t len
    )
{
    uint32_t aligned_len = (len + 0xF) & (~0xF);
    stackinstant.stack_current_alloc_size -= aligned_len;
    stackinstant.stack_current_address -= aligned_len;
}

uint32_t get_stack_current_alloc_size
    (
    void
    )
{
    return stackinstant.stack_current_alloc_size;
}

void reset_stack_ptr_to_assigned_position
    (
    uint32_t assigned_size
    )
{
    stackinstant.stack_current_address = (int8_t*)stackinstant.stack_starting_address + assigned_size;
    stackinstant.stack_current_alloc_size = assigned_size;
}

double** alloc_matrix
    (
    uint32_t nrow,
    uint32_t ncol
    )
{
    //1-index numbering
    double** m = (double**)alloc_from_stack((nrow + 1) * sizeof(double*));
    m[1] = (double*)alloc_from_stack((nrow * ncol + 1) * sizeof(double));
    uint32_t i = 0;
    for (i = 2; i <= nrow; ++i)
        {
        m[i] = m[i - 1] + ncol;
        }
    return m;
}

void free_matrix
    (
    uint32_t nrow,
    uint32_t ncol
    )
{
    partial_free_from_stack((nrow * ncol + 1) * sizeof(double));
    partial_free_from_stack((nrow + 1) * sizeof(double*));
}

void print_matrix
    (
    double** M,
    const uint32_t startrow,
    const uint32_t endrow,
    const uint32_t startcol,
    const uint32_t endcol,
    const double scale
    )
{
    uint32_t i = 1;
    uint32_t j = 1;
    for (j = startrow; j <= endrow; ++j)
        {
        for (i = startcol; i <= endcol; ++i)
            {
            printf("%.3f ", scale * M[j][i]);
            }
        printf("\n");
        }
    printf("========================================\n");
}

void transpose_matrix
    (
    double** src,
    double** dst,
    uint32_t src_rows,
    uint32_t src_cols
    )
{
    uint32_t i = 1;
    uint32_t j = 1;
    for (j = 1; j <= src_rows; ++j)
        {
        for (i = 1; i <= src_cols; ++i)
            {
            dst[i][j] = src[j][i];
            }
        }
}

void multiply_matrix
    (
    double** a,
    double** b,
    double** c,
    uint32_t m,
    uint32_t n,
    uint32_t l
    )
{
    //c = a*b,   a is m by n, b is n by l, hence c is m by l
    uint32_t i = 1;
    uint32_t j = 1;
    uint32_t k = 1;
    for (i = 1; i <= m; ++i)
        {
        for (j = 1; j <= l; ++j)
            {
            c[i][j] = 0.0;
            for (k = 1; k <= n; ++k)
                {
                c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
}

double* alloc_vector
    (
    uint32_t len
    )
{
    //1-index numbering
    double* v = (double*)alloc_from_stack((len + 1)*sizeof(double));
    return v;
}

void free_vector
    (
    uint32_t len
    )
{
    partial_free_from_stack((len + 1)*sizeof(double));
}

void print_vector
    (
    double* V,
    uint32_t len
    )
{
    uint32_t i = 1;
    for (i = 1; i <= len; ++i)
        {
        printf("%f\n", V[i]);
        }
}

void print_image_float
    (
    const float* src,
    const int32_t SRCW,
    const int32_t wstart,
    const int32_t wend,
    const int32_t hstart,
    const int32_t hend,
    const float scale
    )
{
    int32_t i = 0;
    int32_t j = 0;
    for (j = hstart; j < hend; ++j)
        {
        for (i = wstart; i < wend; ++i)
            {
            printf("%.5f ", src[j * SRCW + i] * scale);
            }
        printf("\n");
        }
    printf("========================================\n");
}