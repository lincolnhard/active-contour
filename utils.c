#include "utils.h"

/* will alloc memory inside */
float* create_gaussian_based_first_derivatives
    (
    float sigma,
    DIRECTION direc,
    int32_t* radius
    )
{
    const int32_t KERNEL_RADIUS = (int32_t)(ceil(fabs(sigma) * 3));
    const int32_t KERNEL_WIDTH = (KERNEL_RADIUS << 1) + 1;
    float* kernel = (float*)alloc_from_stack(KERNEL_WIDTH * KERNEL_WIDTH * sizeof(float));
    float* kerneldata = kernel;
    *radius = KERNEL_RADIUS;
    int32_t xidx = 0;
    int32_t yidx = 0;
    const float TWO_SIGMA_SQR = 2 * sigma * sigma;
    const float TWO_PI_SIGMA_DBLSQR = TWO_SIGMA_SQR * (float)M_PI * sigma * sigma;
    int32_t* paraptr = NULL;
    switch (direc)
        {
        case X:
            paraptr = &xidx;
            break;
        case Y:
        default:
            paraptr = &yidx;
            break;
        }
    for (xidx = KERNEL_RADIUS; xidx >= -KERNEL_RADIUS; --xidx)
        {
        for (yidx = KERNEL_RADIUS; yidx >= -KERNEL_RADIUS; --yidx)
            {
            double multiplier = -exp(-(xidx * xidx + yidx * yidx) / TWO_SIGMA_SQR) / TWO_PI_SIGMA_DBLSQR;
            *kerneldata = (float)(*paraptr * multiplier);
            ++kerneldata;
            }
        }
    return kernel;
}

/* will alloc memory inside */
float* create_gaussian_based_second_derivatives
    (
    float sigma,
    DIRECTION direc,
    int32_t* radius
    )
{
    const int32_t KERNEL_RADIUS = (int32_t)(ceil(fabs(sigma) * 3));
    const int32_t KERNEL_WIDTH = (KERNEL_RADIUS << 1) + 1;
    float* kernel = (float*)alloc_from_stack(KERNEL_WIDTH * KERNEL_WIDTH * sizeof(float));
    float* kerneldata = kernel;
    *radius = KERNEL_RADIUS;
    int32_t xidx = 0;
    int32_t yidx = 0;
    const float TWO_SIGMA_SQR = 2 * sigma * sigma;
    const float TWO_PI_SIGMA_DBLSQR = TWO_SIGMA_SQR * (float)M_PI * sigma * sigma;
    int32_t* paraptr = NULL;
    switch (direc)
        {
        case X:
            paraptr = &xidx;
            break;
        case Y:
        default:
            paraptr = &yidx;
            break;
        }
    for (xidx = KERNEL_RADIUS; xidx >= -KERNEL_RADIUS; --xidx)
        {
        for (yidx = KERNEL_RADIUS; yidx >= -KERNEL_RADIUS; --yidx)
            {
            double multiplier = exp(-(xidx * xidx + yidx * yidx) / TWO_SIGMA_SQR) / TWO_PI_SIGMA_DBLSQR;
            *kerneldata = (float)((((*paraptr) * (*paraptr) / (sigma * sigma)) - 1 ) * multiplier);
            ++kerneldata;
            }
        }
    return kernel;
}

void kernel_filter32f
    (
    const float* src,
    const float* kernel,
    const SRCW,
    const SRCH,
    const KERNELR,
    float* dst
    )
{
    int32_t xidx = 0;
    int32_t yidx = 0;
    int32_t idx = 0;
    int32_t kxidx = 0;
    int32_t kyidx = 0;
    float* dstptr = dst;
    for (yidx = 0; yidx < SRCH; ++yidx)
        {
        for (xidx = 0; xidx < SRCW; ++xidx)
            {
            float combinationsum = 0.0f;
            int32_t kidx = 0;
            for (kyidx = -KERNELR; kyidx <= KERNELR; ++kyidx)
                {
                for (kxidx = -KERNELR; kxidx <= KERNELR; ++kxidx)
                    {
                    int32_t tempy = yidx + kyidx;
                    int32_t tempx = xidx + kxidx;
                    if (tempx >= 0 && tempx < SRCW && tempy >= 0 && tempy < SRCH)
                        {
                        int32_t idx = tempx + tempy * SRCW;
                        combinationsum += kernel[kidx] * src[idx];
                        }
                    ++kidx;
                    }
                }
            *dstptr = combinationsum;
            ++dstptr;
            }
        }
}

static double pythag
    (
    double a,
    double b
    )
{
    double absa = fabs(a);
    double absb = fabs(b);
    if (absa > absb)
        {
        return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
        }
    else
        {
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
        }
}

#define SVDSIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x))
static void svdcmp
    (
    double** a,
    uint32_t m,
    uint32_t n,
    double w[],
    double** v
    )
{
    bool flag = 0;
    uint32_t i = 1;
    int its = 1;
    uint32_t j = 1;
    uint32_t jj = 1;
    uint32_t k = 1;
    uint32_t l = 1;
    uint32_t nm = 1;
    double anorm = 0.0;
    double scale = 0.0;
    double g = 0.0;
    double c = 0.0;
    double f = 0.0;
    double h = 0.0;
    double s = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double* rv1 = alloc_vector(n);
    for (i = 1; i <= n; ++i)
        {
        l = i + 1;
        rv1[i] = scale*g;
        g = s = scale = 0.0;
        if (i <= m)
            {
            for (k = i; k <= m; ++k)
                {
                scale += fabs(a[k][i]);
                }
            if (scale)
                {
                for (k = i; k <= m; ++k)
                    {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                    }
                f = a[i][i];
                g = -SVDSIGN(sqrt(s), f);
                h = f*g - s;
                a[i][i] = f - g;
                for (j = l; j <= n; ++j)
                    {
                    for (s = 0.0, k = i; k <= m; ++k)
                        {
                        s += a[k][i] * a[k][j];
                        }
                    f = s / h;
                    for (k = i; k <= m; ++k)
                        {
                        a[k][j] += f*a[k][i];
                        }
                    }
                for (k = i; k <= m; ++k)
                    {
                    a[k][i] *= scale;
                    }
                }
            }
        w[i] = scale *g;
        g = s = scale = 0.0;
        if (i <= m && i != n)
            {
            for (k = l; k <= n; ++k)
                {
                scale += fabs(a[i][k]);
                }
            if (scale)
                {
                for (k = l; k <= n; ++k)
                    {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                    }
                f = a[i][l];
                g = -SVDSIGN(sqrt(s), f);
                h = f*g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; ++k)
                    {
                    rv1[k] = a[i][k] / h;
                    }
                for (j = l; j <= m; ++j)
                    {
                    for (s = 0.0, k = l; k <= n; ++k)
                        {
                        s += a[j][k] * a[i][k];
                        }
                    for (k = l; k <= n; ++k)
                        {
                        a[j][k] += s*rv1[k];
                        }
                    }
                for (k = l; k <= n; ++k)
                    {
                    a[i][k] *= scale;
                    }
                }
            }
        anorm = max(anorm, (fabs(w[i]) + fabs(rv1[i])));
        }
    for (i = n; i >= 1; --i)
        {
        if (i < n)
            {
            if (g)
                {
                for (j = l; j <= n; ++j)
                    {
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                    }
                for (j = l; j <= n; ++j)
                    {
                    for (s = 0.0, k = l; k <= n; ++k)
                        {
                        s += a[i][k] * v[k][j];
                        }
                    for (k = l; k <= n; ++k)
                        {
                        v[k][j] += s*v[k][i];
                        }
                    }
                }
            for (j = l; j <= n; ++j)
                {
                v[i][j] = v[j][i] = 0.0;
                }
            }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
        }
    for (i = min(m, n); i >= 1; --i)
        {
        l = i + 1;
        g = w[i];
        for (j = l; j <= n; ++j)
            {
            a[i][j] = 0.0;
            }
        if (g)
            {
            g = 1.0 / g;
            for (j = l; j <= n; ++j)
                {
                for (s = 0.0, k = l; k <= m; ++k)
                    {
                    s += a[k][i] * a[k][j];
                    }
                f = (s / a[i][i]) * g;
                for (k = i; k <= m; ++k)
                    {
                    a[k][j] += f*a[k][i];
                    }
                }
            for (j = i; j <= m; ++j)
                {
                a[j][i] *= g;
                }
            }
        else
            {
            for (j = i; j <= m; ++j)
                {
                a[j][i] = 0.0;
                }
            }
        ++a[i][i];
        }
    for (k = n; k >= 1; --k)
        {
        for (its = 1; its <= 30; ++its)
            {
            flag = 1;
            for (l = k; l >= 1; --l)
                {
                nm = l - 1;
                if ((double)(fabs(rv1[l]) + anorm) == anorm)
                    {
                    flag = 0;
                    break;
                    }
                if ((double)(fabs(w[nm]) + anorm) == anorm)
                    {
                    break;
                    }
                }
            if (flag)
                {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; ++i)
                    {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];
                    if ((double)(fabs(f) + anorm) == anorm)
                        {
                        break;
                        }
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g*h;
                    s = -f*h;
                    for (j = 1; j <= m; ++j)
                        {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y*c + z*s;
                        a[j][i] = z*c - y*s;
                        }
                    }
                }
            z = w[k];
            if (l == k)
                {
                if (z < 0.0)
                    {
                    w[k] = -z;
                    for (j = 1; j <= n; ++j)
                        {
                        v[j][k] = -v[j][k];
                        }
                    }
                break;
                }
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SVDSIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for (j = l; j <= nm; ++j)
                {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s*g;
                g = c*g;

                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; ++jj)
                    {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x*c + z*s;
                    v[jj][i] = z*c - x*s;
                    }
                z = pythag(f, h);
                w[j] = z;
                if (z)
                    {
                    z = 1.0 / z;
                    c = f*z;
                    s = h*z;
                    }
                f = c*g + s*y;
                x = c*y - s*g;
                for (jj = 1; jj <= m; ++jj)
                    {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y*c + z*s;
                    a[jj][i] = z*c - y*s;
                    }
                }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
            }
        }
    free_vector(n);
}

static void svd_to_inverse
    (
    double** U,
    double* w,
    double** V,
    double** inv,
    int32_t M,
    int32_t N
    )
{
    //U is M by N; w is N ; V is N by N; inv is N by M dimensional
    int32_t i = 1;
    int32_t j = 1;
    int32_t k = 1;
    double** tmp_NN = alloc_matrix(N, N);
    double** tmp_w = alloc_matrix(N, N);

    /*tmp_NN = V * 1/w */
    for (i = 1; i <= N; ++i)
        {
        for (j = 1; j <= N; ++j)
            {
            tmp_w[i][j] = 0.0;
            }
        }
    for (i = 1; i <= N; ++i)
        {
        if (w[i] != 0.0)
            {
            tmp_w[i][i] = 1 / w[i];
            }
        else
            {
            tmp_w[i][i] = 0.0;
            }
        }

    multiply_matrix(V, tmp_w, tmp_NN, N, N, N);

    /*inv = tmp_NN * U' */
    for (i = 1; i <= N; ++i)
        {
        for (j = 1; j <= M; ++j)
            {
            inv[i][j] = 0.0;
            for (k = 1; k <= N; ++k)
                {
                inv[i][j] += tmp_NN[i][k] * U[j][k];
                }
            }
        }
    free_matrix(N, N);
    free_matrix(N, N);
}

void pinv
    (
    double* const* A,
    double** invA,
    uint32_t M,
    uint32_t N
    )
{
    //1-index numbering
    //invA is inverse of A, A is M by N, invA is N by M
    double wmin = 0.0;
    double wmax = 0.0;
    uint32_t i = 1;
    uint32_t j = 1;
    double** U = alloc_matrix(M, N);
    double** V = alloc_matrix(M, N);
    double* w = alloc_vector(N);
    //copy A --> U, svdcmp() will modify input matrix
    memcpy(U[1], A[1], (M * N + 1) * sizeof(double));
    svdcmp(U, M, N, w, V);
    for (i = 1; i <= N; ++i)
        {
        if (w[i] > wmax)
            {
            wmax = w[i];
            }
        }
    const double PINV_SMALL_NUMBER = 1.0e-16;
    wmin = wmax * PINV_SMALL_NUMBER;
    for (i = 1; i <= N; ++i)
        {
        if (w[i] < wmin)
            {
            w[i] = 0.0;
            }
        }
    svd_to_inverse(U, w, V, invA, M, N);
    free_vector(N);
    free_matrix(M, N);
    free_matrix(M, N);
}

void solve_linear
    (
    double** A,
    double** B,
    double** dst,
    uint32_t m,
    uint32_t n
    )
{
    //A: input matrix on the left-hand side of the system
    //B: input matrix on the right-hand side of the system
    //dst: output solution
    //m: rows of A
    //n: cols of A
    double** A_prime = alloc_matrix(n, m); //A'
    double** A_gram = alloc_matrix(n, n); //A'A
    double** A_gram_inv = alloc_matrix(n, n); //inv(A'A)
    double** A_gram_inv_A_prime = alloc_matrix(n, m); //inv(A'A)A'
    transpose_matrix(A, A_prime, m, n);
    multiply_matrix(A_prime, A, A_gram, n, m, n);
    pinv(A_gram, A_gram_inv, n, n);
    multiply_matrix(A_gram_inv, A_prime, A_gram_inv_A_prime, n, n, m);
    multiply_matrix(A_gram_inv_A_prime, B, dst, n, m, 1);
    free_matrix(n, m);
    free_matrix(n, n);
    free_matrix(n, n);
    free_matrix(n, m);
}