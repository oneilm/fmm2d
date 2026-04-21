#pragma once

/*
 * Plan-based Laplace FMM interface (cfmm2d, lfmm2d, rfmm2d).
 *
 * When source/target geometry is fixed and charges vary between calls,
 * build a plan once and reuse it to skip the tree build and interaction-list
 * computation on every subsequent call.
 *
 * At most one plan exists at a time (shared by cfmm2d, lfmm2d, rfmm2d).
 * cfmm2d_destroy_plan frees the plan regardless of which build routine
 * created it.
 *
 * Fortran calling convention: all scalars are passed by pointer, arrays by
 * pointer to the first element.  On Linux/macOS with gfortran the mangled
 * name of subroutine FOO is foo_ (lower-case, trailing underscore).
 *
 * Example (C):
 *   double eps = 1e-9;
 *   int ns = ..., nt = ..., ier;
 *   cfmm2d_build_plan_(&eps, &ns, sources, &nt, targ, &ier);
 *
 *   int nd = 1, ifcharge = 1, ifdipole = 0, iper = 0;
 *   int ifpgh = 1, ifpghtarg = 0;
 *   cfmm2d_execute_plan_(&nd, &ifcharge, charge, &ifdipole, dipstr,
 *       &iper, &ifpgh, pot, grad, hess,
 *       &ifpghtarg, pottarg, gradtarg, hesstarg, &ier);
 *
 *   cfmm2d_destroy_plan_(&ier);
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdint.h>

/* ------------------------------------------------------------------
 * cfmm2d plan  (Cauchy FMM, complex-valued charges / dipoles)
 * ------------------------------------------------------------------ */

/* Build the plan for a given source/target geometry. */
void cfmm2d_build_plan_(
    double       *eps,
    int          *ns,
    double       *sources,   /* (2,ns) */
    int          *nt,
    double       *targ,      /* (2,nt) */
    int          *ier);

/*
 * Execute the FMM with the pre-built plan.
 * charge(nd,ns), dipstr(nd,ns) : complex source densities
 * pot/grad/hess(nd,ns)         : Cauchy potential/gradient/hessian at sources
 * pottarg/gradtarg/hesstarg    : same at targets
 */
void cfmm2d_execute_plan_(
    int          *nd,
    int          *ifcharge,
    double _Complex *charge,    /* (nd,ns) */
    int          *ifdipole,
    double _Complex *dipstr,    /* (nd,ns) */
    int          *iper,
    int          *ifpgh,
    double _Complex *pot,       /* (nd,ns) */
    double _Complex *grad,      /* (nd,ns) */
    double _Complex *hess,      /* (nd,ns) */
    int          *ifpghtarg,
    double _Complex *pottarg,   /* (nd,nt) */
    double _Complex *gradtarg,  /* (nd,nt) */
    double _Complex *hesstarg,  /* (nd,nt) */
    int          *ier);

/* Free all plan memory. */
void cfmm2d_destroy_plan_(int *ier);


/* ------------------------------------------------------------------
 * lfmm2d plan  (complex-charge Laplace FMM)
 * ------------------------------------------------------------------ */

/* Same geometry as cfmm2d; delegates to cfmm2d_build_plan. */
void lfmm2d_build_plan_(
    double       *eps,
    int          *ns,
    double       *sources,   /* (2,ns) */
    int          *nt,
    double       *targ,      /* (2,nt) */
    int          *ier);

/*
 * charge(nd,ns), dipstr(nd,ns) : complex densities
 * dipvec(nd,2,ns)              : real dipole orientation vectors
 * pot(nd,ns)                   : Laplace potential at sources
 * grad(nd,2,ns)                : gradient [d/dx, d/dy] at sources
 * hess(nd,3,ns)                : hessian [d^2/dx^2, d^2/dxdy, d^2/dy^2]
 */
void lfmm2d_execute_plan_(
    int          *nd,
    int          *ifcharge,
    double _Complex *charge,    /* (nd,ns) */
    int          *ifdipole,
    double _Complex *dipstr,    /* (nd,ns) */
    double       *dipvec,       /* (nd,2,ns) */
    int          *iper,
    int          *ifpgh,
    double _Complex *pot,       /* (nd,ns)   */
    double _Complex *grad,      /* (nd,2,ns) */
    double _Complex *hess,      /* (nd,3,ns) */
    int          *ifpghtarg,
    double _Complex *pottarg,   /* (nd,nt)   */
    double _Complex *gradtarg,  /* (nd,2,nt) */
    double _Complex *hesstarg,  /* (nd,3,nt) */
    int          *ier);


/* ------------------------------------------------------------------
 * rfmm2d plan  (real-charge Laplace FMM)
 * ------------------------------------------------------------------ */

/* Same geometry as cfmm2d; delegates to cfmm2d_build_plan. */
void rfmm2d_build_plan_(
    double       *eps,
    int          *ns,
    double       *sources,   /* (2,ns) */
    int          *nt,
    double       *targ,      /* (2,nt) */
    int          *ier);

/*
 * charge(nd,ns), dipstr(nd,ns) : real densities
 * dipvec(nd,2,ns)              : real dipole orientation vectors
 * pot(nd,ns)                   : Laplace potential at sources
 * grad(nd,2,ns)                : gradient at sources
 * hess(nd,3,ns)                : hessian at sources
 */
void rfmm2d_execute_plan_(
    int    *nd,
    int    *ifcharge,
    double *charge,    /* (nd,ns) */
    int    *ifdipole,
    double *dipstr,    /* (nd,ns) */
    double *dipvec,    /* (nd,2,ns) */
    int    *iper,
    int    *ifpgh,
    double *pot,       /* (nd,ns)   */
    double *grad,      /* (nd,2,ns) */
    double *hess,      /* (nd,3,ns) */
    int    *ifpghtarg,
    double *pottarg,   /* (nd,nt)   */
    double *gradtarg,  /* (nd,2,nt) */
    double *hesstarg,  /* (nd,3,nt) */
    int    *ier);

/* cfmm2d_destroy_plan_ also frees lfmm2d and rfmm2d plans. */

#ifdef __cplusplus
}
#endif
