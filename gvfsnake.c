#include "gvfsnake.h"

void calc_external_force
    (
    const float* src,
    float* u,           /* dEext/dx */
    float* v,           /* dEext/dy */
    const int32_t SRCWIDTH,
    const int32_t SRCHEIGHT
    )
{
    uint32_t current_allocsize = get_stack_current_alloc_size();
    const int32_t SRCSIZE = SRCWIDTH * SRCHEIGHT;
    const int32_t SRCBYTE = SRCSIZE * sizeof(float);
    int32_t gdx_radius = 0;
    float* GDX = create_gaussian_based_first_derivatives(SIGMA_ONE, X, &gdx_radius);
    int32_t gdy_radius = 0;
    float* GDY = create_gaussian_based_first_derivatives(SIGMA_ONE, Y, &gdy_radius);
    int32_t gdxx_radius = 0;
    float* GDXX = create_gaussian_based_second_derivatives(SIGMA_TWO, X, &gdxx_radius);
    int32_t gdyy_radius = 0;
    float* GDYY = create_gaussian_based_second_derivatives(SIGMA_TWO, Y, &gdyy_radius);
    float* buf1 = (float*)alloc_from_stack(SRCBYTE);
    float* buf2 = (float*)alloc_from_stack(SRCBYTE);
    float* buf3 = (float*)alloc_from_stack(SRCBYTE);
    kernel_filter32f(src, GDX, SRCWIDTH, SRCHEIGHT, gdx_radius, buf1); /* I->Ix */
    kernel_filter32f(src, GDY, SRCWIDTH, SRCHEIGHT, gdy_radius, buf2); /* I->Iy */
    int32_t i = 0;
    int32_t j = 0;
    for (i = 0; i < SRCSIZE; ++i)
        {
        buf3[i] = (float)(-2.0 * sqrt(buf1[i] * buf1[i] + buf2[i] * buf2[i])); /* E external aka. f */
        }
    kernel_filter32f(buf3, GDX, SRCWIDTH, SRCHEIGHT, gdx_radius, buf1); /* f->fx */
    kernel_filter32f(buf3, GDY, SRCWIDTH, SRCHEIGHT, gdy_radius, buf2); /* f->fy */
    const float TWO_SIGMA_SQR = 2 * SIGMA_ONE * SIGMA_ONE;
    for (i = 0; i < SRCSIZE; ++i)
        {
        buf1[i] *= -TWO_SIGMA_SQR; /* fx */
        buf2[i] *= -TWO_SIGMA_SQR; /* fy */
        }
    /* calculate GVF (gradient vector flow) */
    for (i = 0; i < SRCSIZE; ++i)
        {
        buf3[i] = buf1[i] * buf1[i] + buf2[i] * buf2[i]; /* magnitude */
        }
    memcpy(u, buf1, SRCBYTE);
    memcpy(v, buf2, SRCBYTE);
    float* uLoG = (float*)alloc_from_stack(SRCBYTE);
    float* vLoG = (float*)alloc_from_stack(SRCBYTE);
    float* GDXXYY = GDXX;
    const int32_t LOG_KERNELSIZE = ((gdxx_radius << 1) + 1) * (gdxx_radius << 1) + 1;
    for (i = 0; i < LOG_KERNELSIZE; ++i)
        {
        GDXXYY[i] = GDXX[i] + GDYY[i]; /* LoG kernel */
        }
    for (j = 0; j < ITER_TIMES_EXT; ++j)
        {
        kernel_filter32f(u, GDXXYY, SRCWIDTH, SRCHEIGHT, gdxx_radius, uLoG);
        kernel_filter32f(v, GDXXYY, SRCWIDTH, SRCHEIGHT, gdxx_radius, vLoG);
        for (i = 0; i < SRCSIZE; ++i)
            {
            u[i] = u[i] + MU * uLoG[i] - buf3[i] * (u[i] - buf1[i]);
            v[i] = v[i] + MU * vLoG[i] - buf3[i] * (v[i] - buf2[i]);
            }
        }
    reset_stack_ptr_to_assigned_position(current_allocsize);
}

double** create_internal_force_pentadiagonal_matrix
    (
    const int32_t numpts
    )
{
    double c0 = 2 * ALPHA + 6 * BETA;
    double c1 = -ALPHA - 4 * BETA;
    double c2 = BETA;
    double** A_inv = alloc_matrix(numpts, numpts);
    double** A = alloc_matrix(numpts, numpts);
    int32_t i = 0;
    for (i = 1; i <= numpts; ++i)
        {
        A[i][i] = c0 + GAMMA;
        }
    for (i = 1; i < numpts; ++i)
        {
        A[i + 1][i] = c1;
        A[i][i + 1] = c1;
        }
    for (i = 1; i < numpts - 1; ++i)
        {
        A[i + 2][i] = c2;
        A[i][i + 2] = c2;
        }
    A[numpts][1] = A[1][numpts] = c1;
    A[numpts - 1][1] = A[numpts][2] = A[1][numpts - 1] = A[2][numpts] = c2;
    pinv(A, A_inv, numpts, numpts);
    free_matrix(numpts, numpts);
    return A_inv;
}

void contour_update
    (
    double** ptsx,
    double** ptsy,
    float* Ex,
    float* Ey,
    double** Eint_coefficients,
    const int32_t numpts,
    const int32_t IMW,
    const int32_t IMH
    )
{
    double** ssx = alloc_matrix(numpts, 1);
    double** ssy = alloc_matrix(numpts, 1);
    int32_t i = 0;
    for (i = 1; i <= numpts; ++i)
        {
        int32_t idx = IMW * (int32_t)(ptsy[i][1]) + (int32_t)(ptsx[i][1]);
        ssx[i][1] = GAMMA * ptsx[i][1] + 4 * Ey[idx];
        ssy[i][1] = GAMMA * ptsy[i][1] + 4 * Ex[idx];
        }
    multiply_matrix(Eint_coefficients, ssx, ptsx, numpts, numpts, 1);
    multiply_matrix(Eint_coefficients, ssy, ptsy, numpts, numpts, 1);
    free_matrix(numpts, 1);
    free_matrix(numpts, 1);
}
