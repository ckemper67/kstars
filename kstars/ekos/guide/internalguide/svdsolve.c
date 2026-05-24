/*
**  - - - - - - - - - - -
**   s v d s o l v e . c
**  - - - - - - - - - - -
**
**  SVD-based linear least-squares solver.
**
**  Implements the Golub-Reinsch algorithm (bidiagonalization + QR
**  iteration with implicit Wilkinson shift) for computing the thin SVD
**  of a tall m x n matrix (m >= n), then uses it to solve overdetermined
**  systems in the least-squares sense.
**
**  References:
**    Golub, G.H., Reinsch, C., 1970, "Singular value decomposition and
**    least squares solutions", Numer. Math. 14, 403-420.
**
**    Golub, G.H., Van Loan, C.F., 2013, Matrix Computations, 4th ed.,
**    Johns Hopkins University Press.
*/

#define SVD_MAX_COLS 64
#include <math.h>

/*
** Helper: 2-norm of (a, b) without undue overflow/underflow.
*/
static double svd_hypot(double a, double b)
{
    double r;
    a = fabs(a);
    b = fabs(b);
    if ( a > b ) {
        r = b / a;
        return a * sqrt(1.0 + r*r);
    }
    if ( b == 0.0 ) return 0.0;
    r = a / b;
    return b * sqrt(1.0 + r*r);
}

/*
** - - - - - - - - -
**  s v d
** - - - - - - - - -
**
**  Thin SVD of an m x n matrix A (m >= n): A = U * diag(w) * V^T.
**
**  Given:
**     m    int       number of rows (observations * 2)
**     n    int       number of columns (model terms)
**     a    double*   m x n matrix, row-major (OVERWRITTEN with U on exit)
**     w    double*   output: n singular values (non-negative)
**     v    double*   output: n x n matrix V (NOT V^T)
**
**  Returned (function value):
**           int      0  = success
**            -1          = illegal dimensions (m < n or n < 1)
**                          singular value (1-based index). Results are
**                          usually trustworthy; treat as a warning and
**                          proceed normally.
**
**  Notes:
**  1) On entry a[i*n+j] is element (i,j).  On exit a contains U, with
**     U^T U = I (n x n identity).
**  2) Singular values in w[] are non-negative but NOT sorted.
**  3) v[i*n+j] is element (i,j) of V, so the full decomposition is
**     a_original = U * diag(w) * v^T.
*/
int svd(int m, int n, double *a, double *w, double *v)
{
    double rv1[SVD_MAX_COLS];
    int i, j, k, l, its, jj, nm;
    int flag, jstat;
    double anorm, c, f, g, h, s, scale, x, y, z;

    if ( n < 1 || m < n || n > SVD_MAX_COLS ) return -1;
    jstat = 0;

    /* ---------------------------------------------------------------
    ** Phase 1: Householder bidiagonalization.
    ** Reduces a to upper bidiagonal form stored in diag(w) and
    ** superdiag(rv1).  Transformations accumulate in a (left) and v (right).
    ** --------------------------------------------------------------- */

    g = 0.0; scale = 0.0; anorm = 0.0;

    for ( i = 0; i < n; i++ ) {
        l = i + 1;
        rv1[i] = scale * g;
        g = 0.0; s = 0.0; scale = 0.0;

        /* Left Householder: zero below diagonal in column i */
        for ( k = i; k < m; k++ ) scale += fabs(a[k*n+i]);
        if ( scale != 0.0 ) {
            for ( k = i; k < m; k++ ) {
                a[k*n+i] /= scale;
                s += a[k*n+i] * a[k*n+i];
            }
            f = a[i*n+i];
            g = -copysign(sqrt(s), f);
            h = f*g - s;
            a[i*n+i] = f - g;
            for ( j = l; j < n; j++ ) {
                s = 0.0;
                for ( k = i; k < m; k++ ) s += a[k*n+i] * a[k*n+j];
                f = s / h;
                for ( k = i; k < m; k++ ) a[k*n+j] += f * a[k*n+i];
            }
            for ( k = i; k < m; k++ ) a[k*n+i] *= scale;
        }
        w[i] = scale * g;

        g = 0.0; s = 0.0; scale = 0.0;

        /* Right Householder: zero right of superdiagonal in row i */
        if ( i < m && i != n-1 ) {
            for ( k = l; k < n; k++ ) scale += fabs(a[i*n+k]);
            if ( scale != 0.0 ) {
                for ( k = l; k < n; k++ ) {
                    a[i*n+k] /= scale;
                    s += a[i*n+k] * a[i*n+k];
                }
                f = a[i*n+l];
                g = -copysign(sqrt(s), f);
                h = f*g - s;
                a[i*n+l] = f - g;
                for ( k = l; k < n; k++ ) rv1[k] = a[i*n+k] / h;
                for ( j = l; j < m; j++ ) {
                    s = 0.0;
                    for ( k = l; k < n; k++ ) s += a[j*n+k] * a[i*n+k];
                    for ( k = l; k < n; k++ ) a[j*n+k] += s * rv1[k];
                }
                for ( k = l; k < n; k++ ) a[i*n+k] *= scale;
            }
        }
        anorm = fmax(anorm, fabs(w[i]) + fabs(rv1[i]));
    }

    /* ---------------------------------------------------------------
    ** Phase 2: Accumulate right-hand transformations into V.
    ** --------------------------------------------------------------- */

    g = 0.0;
    for ( i = n-1; i >= 0; i-- ) {
        l = i + 1;
        if ( i < n-1 ) {
            if ( g != 0.0 ) {
                for ( j = l; j < n; j++ )
                    v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
                for ( j = l; j < n; j++ ) {
                    s = 0.0;
                    for ( k = l; k < n; k++ ) s += a[i*n+k] * v[k*n+j];
                    for ( k = l; k < n; k++ ) v[k*n+j] += s * v[k*n+i];
                }
            }
            for ( j = l; j < n; j++ ) v[i*n+j] = v[j*n+i] = 0.0;
        }
        v[i*n+i] = 1.0;
        g = rv1[i];
    }

    /* ---------------------------------------------------------------
    ** Phase 3: Accumulate left-hand transformations into U (in-place in a).
    ** --------------------------------------------------------------- */

    for ( i = n-1; i >= 0; i-- ) {
        l = i + 1;
        g = w[i];
        for ( j = l; j < n; j++ ) a[i*n+j] = 0.0;
        if ( g != 0.0 ) {
            g = 1.0 / g;
            for ( j = l; j < n; j++ ) {
                s = 0.0;
                for ( k = l; k < m; k++ ) s += a[k*n+i] * a[k*n+j];
                f = (s / a[i*n+i]) * g;
                for ( k = i; k < m; k++ ) a[k*n+j] += f * a[k*n+i];
            }
            for ( j = i; j < m; j++ ) a[j*n+i] *= g;
        } else {
            for ( j = i; j < m; j++ ) a[j*n+i] = 0.0;
        }
        a[i*n+i] += 1.0;
    }

    /* ---------------------------------------------------------------
    ** Phase 4: Diagonalization of bidiagonal form.
    ** QR iteration with implicit Wilkinson shift, cycling over singular values
    ** from largest index to smallest.
    ** --------------------------------------------------------------- */

    for ( k = n-1; k >= 0; k-- ) {

        for ( its = 0; its < 30; its++ ) {

            /* Test for splitting.  Look for a small superdiagonal element. */
            flag = 1;
            l = k;
            while ( l > 0 ) {
                nm = l - 1;
                if ( fabs(rv1[l]) + anorm == anorm ) { flag = 0; break; }
                if ( fabs(w[nm])  + anorm == anorm ) break;
                l--;
            }

            if ( flag ) {
                /* rv1[l] is nonzero and l > 0: cancel rv1[l] by
                   Givens rotation in the (nm, l) plane. */
                c = 0.0; s = 1.0;
                for ( i = l; i <= k; i++ ) {
                    f = s * rv1[i];
                    rv1[i] *= c;
                    if ( fabs(f) + anorm == anorm ) break;
                    g = w[i];
                    h = svd_hypot(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c =  g * h;
                    s = -f * h;
                    for ( j = 0; j < m; j++ ) {
                        y = a[j*n+nm];
                        z = a[j*n+i];
                        a[j*n+nm] = y*c + z*s;
                        a[j*n+i]  = z*c - y*s;
                    }
                }
            }

            z = w[k];

            if ( l == k ) {
                /* Convergence: ensure singular value is non-negative. */
                if ( z < 0.0 ) {
                    w[k] = -z;
                    for ( j = 0; j < n; j++ ) v[j*n+k] = -v[j*n+k];
                }
                break;
            }

            if ( its == 29 ) {
                /* Did not fully converge: record warning, continue. */
                jstat = k + 1;   /* 1-based index of problematic singular value */
            }

            /* Compute implicit shift from bottom 2x2 minor. */
            x  = w[l];
            nm = k - 1;
            y  = w[nm];
            g  = rv1[nm];
            h  = rv1[k];
            f  = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y);
            g  = svd_hypot(f, 1.0);
            f  = ((x-z)*(x+z) + h*((y / (f + copysign(g, f))) - h)) / x;

            /* Next QR transformation. */
            c = 1.0; s = 1.0;
            for ( j = l; j <= nm; j++ ) {
                i  = j + 1;
                g  = rv1[i];
                y  = w[i];
                h  = s*g;
                g  = c*g;
                z  = svd_hypot(f, h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;
                for ( jj = 0; jj < n; jj++ ) {
                    x = v[jj*n+j];
                    z = v[jj*n+i];
                    v[jj*n+j] = x*c + z*s;
                    v[jj*n+i] = z*c - x*s;
                }
                z = svd_hypot(f, h);
                w[j] = z;
                if ( z != 0.0 ) {
                    z = 1.0 / z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for ( jj = 0; jj < m; jj++ ) {
                    y = a[jj*n+j];
                    z = a[jj*n+i];
                    a[jj*n+j] = y*c + z*s;
                    a[jj*n+i] = z*c - y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k]   = x;
        }
    }

    return jstat;
}

/*
** - - - - - - - - - - - - - - -
**  s v d b a c k s u b
** - - - - - - - - - - - - - - -
**
**  Least-squares solve using pre-computed SVD: A = U * diag(w) * V^T.
**
**  Given:
**     m    int       rows of A (= rows of U)
**     n    int       columns of A (= length of w, size of V)
**     u    double*   m x n matrix U (from svd, row-major)
**     w    double*   n singular values (from svd)
**     v    double*   n x n matrix V (from svd, row-major)
**     b    double*   m-vector (right-hand side)
**     wmin double    threshold: singular values <= wmin are treated as zero
**
**  Returned:
**     x    double*   n-vector (least-squares solution, minimum-norm)
**
**  Note: singular values at or below wmin are zeroed in the inverse,
**  giving the minimum-norm solution when A is rank-deficient.
*/
void svd_backsub(int m, int n,
                      const double *u, const double *w, const double *v,
                      const double *b, double wmin,
                      double *x)
{
    double tmp[SVD_MAX_COLS];
    double s;
    int i, j;

    /* tmp = diag(1/w) * U^T * b, zeroing near-singular components */
    for ( j = 0; j < n; j++ ) {
        s = 0.0;
        if ( w[j] > wmin ) {
            for ( i = 0; i < m; i++ ) s += u[i*n+j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }

    /* x = V * tmp */
    for ( j = 0; j < n; j++ ) {
        s = 0.0;
        for ( i = 0; i < n; i++ ) s += v[j*n+i] * tmp[i];
        x[j] = s;
    }
}
