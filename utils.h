#if !defined UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif
#define _USE_MATH_DEFINES

#include "DS.h"
#include <math.h>

/* will alloc memory inside */
float*  create_gaussian_based_first_derivatives
    (
    float sigma,
    DIRECTION direc,
    int32_t* radius
    );

/* will alloc memory inside */
float* create_gaussian_based_second_derivatives
    (
    float sigma,
    DIRECTION direc,
    int32_t* radius
    );

void kernel_filter32f
    (
    const float* src,
    const float* kernel,
    const int32_t srcw,
    const int32_t srch,
    const int32_t kernelradius,
    float* dst
    );

void pinv
    (
    double* const* A,
    double** invA,
    uint32_t M,
    uint32_t N
    );

void solve_linear
    (
    double** A,
    double** B,
    double** dst,
    uint32_t m,
    uint32_t n
    );

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */