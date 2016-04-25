#ifndef __EOS_H__
#define __EOS_H__

#include <stddef.h>
#include <math.h>
#include "libguile.h"
#include "eos.h"

/*
  Function to calculate omega_ai(T)
*/
typedef double (OmegaAiFunc) (CubicEosModel*, size_t, double);

struct CubicEosModel {
  /*
    +--------------------+--------------------+--------------------+
    |EOS                 |m1                  |m2                  |
    +--------------------+--------------------+--------------------+
    |Redlich-Kwong (RK)  |0                   |1                   |
    +--------------------+--------------------+--------------------+
    |Soave-RK (SRK)      |0                   |1                   |
    +--------------------+--------------------+--------------------+
    |Peng-Robinson (PR)  |1+sqrt(2)           |1-sqrt(2)           |
    +--------------------+--------------------+--------------------+
  */
  SCM type;                     /* The type of EOS, symbol (PR/PRCORR/RK/SRK) */
  double m1;
  double m2;

  size_t nc;                    /* number of components */
  SCM cnames;                   /* component names, vector of symbols */
  double* mw;                   /* molecular weight */
  double* pc;                   /* critical pressure */
  double* tc;                   /* critical temperatureN */
  double* w;                    /* eccentric factor */
  double* bic;                  /* binary interaction coefficient */
  /*
    +--------------------+--------------------+--------------------+
    |EOS                 |omega_a0            |omega_b0            |
    +--------------------+--------------------+--------------------+
    |RK SRK              |0.4274802           |0.08664035          |
    +--------------------+--------------------+--------------------+
    |PR PRCORR           |0.457235529         |0.07796074          |
    +--------------------+--------------------+--------------------+
  */
  double* omega_a0;             /* override default value */
  double* omega_b0;             /* override default value */
  double* vshift;               /* volume shift */

  OmegaAiFunc* omega_ai_func;

  /* -------------------------Working Zone-------------------------------------- */
  double* poly_w;               /* polynomial of w_i in omega_ai calculation */
  size_t rank;                  /* rank used for reduced variable calculation */
  double* eigv;                 /* (rank) eigenvectors of (1-bic) */
  double* lambda;               /* (rank) eigenvalues of (1-bic)*/
  double t;                     /* current temperature */
  double p;                     /* current pressure */
  double* sq_ai_t;              /* sqrt(omega_ai(t))/t */
  double* sq_ai;                /* sqrt(ai) = sq_ai_t * sqrt(pri) */
  double* bi_t;                 /* omega_bi/tri */
  double* bi;                   /* bi = bi_t * pri */
  double* q_mat;                /* (rank+1)*nc expanded transformation matrix */
};

static double omega_ai_func_pr (CubicEosModel* model, size_t i, double t);
static double omega_ai_func_prcorr (CubicEosModel* model, size_t i, double t);
static double omega_ai_func_srk (CubicEosModel* model, size_t i, double t);
static double omega_ai_func_rk (CubicEosModel* model, size_t i, double t);
static void init_cubic_eos_model_pr (CubicEosModel* model);
static void init_cubic_eos_model_prcorr (CubicEosModel* model);
static void init_cubic_eos_model_rk (CubicEosModel* model);
static void init_cubic_eos_model_srk (CubicEosModel* model);
static void init_cubic_eos_model (CubicEosModel* model);
static void update_t(CubicEosModel* model, double t);
static void update_p(CubicEosModel* model, double p);
static void update_tp (CubicEosModel* model, double t, double p);


int dspev (char JOBZ, char UPLO,
           int N, double *A,
           double *W, double *Z,
           int LDZ, double *WORK);


#endif  /* __EOS_H__ */
