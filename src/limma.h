#ifndef LIMMA_H
#define LIMMA_H

typedef struct {
	double *x;      /* compact design matrix -> dqrls QR storage (nobs x nbeta) */
	double *y;      /* compact response vector (nobs) */
	double *coeff;  /* dqrls coefficients in pivot order (nbeta) */
	double *resid;  /* dqrls residual working output (nobs) */
	double *effec;  /* rotated response Q^T y from dqrls (nobs) */
	double *qraux;  /* dqrls Householder auxiliary data (nbeta) */
	double *work;   /* dqrls workspace (2 * nbeta) */
	double *rinv;   /* compact inverse of estimable R factor (rank x rank) */
	double *cest;   /* contrast restricted to estimable columns in pivot order */
	double *scale;  /* per-observation sqrt(weight), used by GLS whitening */
	int *obser;     /* original array index for each compact observation (nobs) */
	int *pivot;     /* 1-based dqrls pivot vector (nbeta) */
	int *estpos;    /* original-column to pivot-position map for contrasts */
} lmws;

typedef struct {
	lmws qr;    /* ordinary least-squares workspace after GLS prewhitening */
	double *v;  /* observed covariance submatrix -> Cholesky factor U */
} glsws;

typedef struct {
	double *x;      /* weighted fixed-effect design -> dqrls QR storage */
	double *rhs;    /* dqrls RHS [block indicators | response], then gamma data */
	double *coeff;  /* dqrls coefficient output for all RHS columns */
	double *resid;  /* dqrls residual working output for all RHS columns */
	double *effec;  /* Q^T [Z | y], including REML residual-space rows */
	double *qraux;  /* dqrls Householder auxiliary data */
	double *qwork;  /* dqrls workspace */
	double *qtz;    /* Q2^T Z -> LAPACK QR storage (mq x nblocks) */
	double *tau;    /* DGEQRF reflector scalars for qtz */
	double *rmat;   /* top triangular QR block copied for small DGESVD */
	double *svals;  /* singular values of Q2^T Z from the small SVD */
	double *u;      /* left singular vectors used to rotate residual components */
	double *gvec;   /* restricted response Q2^T y -> Q^T(Q2^T y) */
	double *swork;  /* LAPACK workspace shared by QR, Q multiply, and SVD */
	double *mval;   /* compact finite expression values for one gene */
	double *scale;  /* compact sqrt(weight) factors for one gene */
	int *obser;     /* original array index for each compact observation */
	int *bmap;      /* observed block-level relabeling to 0..nblocks-1 */
	int *pivot;     /* 1-based fixed-effect dqrls pivot vector */
	int lwork;      /* allocated length of swork */
} dupws;

typedef struct {
	double *xw;     /* weighted design -> QR reflectors -> Q1 (narrays x p) */
	double *yw;     /* weighted response (narrays) */
	double *rw;     /* residual vector (narrays) */
	double *tau;    /* QR reflector scalars (p) */
	double *cvec;   /* Q1^T yw (p) */
	double *q2;     /* column products of Q1 (narrays x p2) */
	double *bmat;   /* Q2^T Z  (p2 x (ngam+1)) */
	double *hvec;   /* leverages diag(Q1 Q1^T) (narrays) */
	double *info;   /* Fisher info incl. intercept ((ngam+1)^2) */
	double *pinfo2; /* per-thread accumulator for info2 (ngam^2) */
	double *pz;     /* per-thread accumulator for z (narrays) */
	double *work;   /* LAPACK workspace */
	int lwork;      /* allocated length of work */
} awws;

typedef struct {
	double *y;      /* count response for one gene (narrays) */
	double *uw;     /* weighted IRLS working response sqrt(mu) * u */
	double *eta;    /* linear predictor offset + design beta (narrays) */
	double *mu;     /* fitted Poisson means exp(eta) (narrays) */
	double *xw;     /* weighted IRLS design sqrt(mu) * design -> QR storage */
	double *coef;   /* dqrls coefficients in pivot order (p) */
	double *beta;   /* unpivoted IRLS coefficients in original column order */
	double *resid;  /* dqrls residual working output (narrays) */
	double *effec;  /* dqrls Q^T uw working output (narrays) */
	double *qraux;  /* dqrls Householder auxiliary data (p) */
	double *work;   /* dqrls workspace (2 * p) */
	int *pivot;     /* 1-based dqrls pivot vector (p) */
} poisws;

typedef struct {
	int *seed;       /* sampling-point indices (sized npts; nseeds <= npts) */
	int *fstart;     /* span start per seed */
	int *fend;       /* span end per seed */
	double *fdist;   /* max covariate distance per seed */
	double *work;    /* per-point weights / |residual| scratch (npts) */
	int *ror;        /* sort-order index (npts) */
} lowessws;

int clampthreads(int, int);

int qrfit(int, int, int, int, const double *, int, double *, double *, double *, int *, lmws *);

int lmgene(const double *, const double *, const double *, int, int, int, const double *, int, int, double *, double *, double *, int *, lmws *);

int glsgene(const double *, const double *, const double *, const double *, int, int, int, const double *, int, int, double *, double *, double *, int *, const double *, glsws *);

double dupcorgene(const double *, const double *, const int *, const double *, int, int, int, int, int, dupws *);

void lm(const double *, const double *, const double *, int, int, int, const double *, int, int, double *, double *, double *, int *);

int gls(const double *, const double *, const double *, const double *, int, int, int, const double *, int, int, double *, double *, double *, int *);

int dupcor(const double *, const double *, const int *, const double *, int, int, int, int, int, double *);

void poisgene(const double *, const double *, const double *, const double *, int, int, int, int, double *, poisws *);

void pois(const double *, const double *, const double *, const double *, int, int, int, int, double *);

void lowess(const double *, const double *, const double *, int, double, int, double, double *, double *);

void awremlgene(const double *, const double *, const double *, const double *, int, int, int, int, int, const double *, int, awws *);

int awreml(const double *, const double *, const double *, const double *, int, int, int, int, double, int, double, int, int, double *, int *);

#endif
