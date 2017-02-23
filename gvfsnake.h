#if !defined GVF_H
#define GVF_H

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SIGMA_ONE (8.0f) /* sigma of first order Gaussian derivatives */
#define SIGMA_TWO (1.0f) /* sigma of second order Gaussian derivatives */
#define ITER_TIMES_EXT (600) /* iteration times for external force calulation */
#define MU (0.2f) /* tradeoff scalar between noise and real edge forces, big mu for noisy image, small mu for distanced initial contour */
#define ALPHA (0.1) /* membrame energy */
#define BETA (0.1) /* thin plate energy */
#define GAMMA (1.0) /* step size (time) */
#define ITER_TIMES_CONTOUR (400) /* iteration times for contour movement */

void calc_external_force
    (
    const float* src,
    float* Eextx,           /* dEext/dx */
    float* Eexty,           /* dEext/dy */
    const int32_t width,
    const int32_t height
    );

/* alloc memory */
double** create_internal_force_pentadiagonal_matrix
    (
    const int32_t numpts
    );

void contour_update
    (
    double** ptsx,
    double** ptsy,
    float* Ex,
    float* Ey,
    double** Eint_coefficients,
    const int32_t numpts,
    const int32_t width,
    const int32_t height
    );

#ifdef __cplusplus
}
#endif

#endif /* GVF_H */