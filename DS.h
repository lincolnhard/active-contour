#if !defined DS_H
#define DS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define STACK_TOTAL_SIZE (0xB00000) //11MB

typedef enum
    {
    X,
    Y
    }DIRECTION;

typedef struct
    {
    int32_t x;
    int32_t y;
    }point2ds32i;

typedef struct
    {
    void* stack_starting_address;
    int8_t* stack_current_address;
    uint32_t stack_current_alloc_size;
    }stack_datastruct;

void init_stack
    (
    void
    );

void free_stack
    (
    void
    );

void reset_stack_ptr_to_initial_position
    (
    void
    );

void* alloc_from_stack
    (
    uint32_t len
    );

void partial_free_from_stack
    (
    uint32_t len
    );

uint32_t get_stack_current_alloc_size
    (
    void
    );

void reset_stack_ptr_to_assigned_position
    (
    uint32_t assigned_size
    );

double** alloc_matrix
    (
    uint32_t nrow,
    uint32_t ncol
    );

void free_matrix
    (
    uint32_t nrow,
    uint32_t ncol
    );

void print_matrix
    (
    double** M,
    const uint32_t startrow,
    const uint32_t endrow,
    const uint32_t startcol,
    const uint32_t endcol,
    const double scale
    );

void transpose_matrix
    (
    double** src,
    double** dst,
    uint32_t src_rows,
    uint32_t src_cols
    );

void multiply_matrix
    (
    double** a,
    double** b,
    double** c,
    uint32_t m,
    uint32_t n,
    uint32_t l
    );

double* alloc_vector
    (
    uint32_t len
    );

void free_vector
    (
    uint32_t len
    );

void print_vector
    (
    double* V,
    uint32_t len
    );

void print_image_float
    (
    const float* src,
    const int32_t SRCW,
    const int32_t wstart,
    const int32_t wend,
    const int32_t hstart,
    const int32_t hend,
    const float scale
    );

#endif /* DS_H */