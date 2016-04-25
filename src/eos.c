#include <math.h>
#include "cblas.h"
#include "libguile.h"
#include "eos_i.h"

#define DEBUG

SCM_GLOBAL_VARIABLE_INIT(eos_cubic_rv_cutoff, "eos-cubic-rv-cutoff", scm_from_double(0.001));

static const double GAS_CONST = 10.73159;

static inline int scm_symbol_is_equal (SCM s1, SCM s2)
{
  return scm_is_true (scm_string_eq (scm_symbol_to_string (s1),
                                     scm_symbol_to_string (s2),
                                     SCM_UNDEFINED, SCM_UNDEFINED,
                                     SCM_UNDEFINED, SCM_UNDEFINED));
}

/* Calculate poly_w for different EOS models */
static inline double compute_poly_w_pr (double w)
{
  return 0.37464 + 1.54226*w - 0.26992*w*w;
}

static inline double compute_poly_w_srk (double w)
{
  return 0.48 + 1.574*w - 0.176*w*w;
}

static inline double compute_poly_w_prcorr (double w)
{
  return w > 0.49 ?
    0.379642 + 1.48503*w - 0.164423*w*w + 0.016666*w*w*w :
    0.37464 + 1.54226*w - 0.26992*w*w;
}

/* OmegaAiFunc for different EOS models */
static double omega_ai_func_pr (CubicEosModel* model, size_t i, double t)
{
  double sqr_tri = sqrt(t/model->tc[i]);
  double fw = model->poly_w[i];
  double oma0 = model->omega_a0[i];
  double c = 1.0 + fw * (1.0 - sqr_tri);
  return oma0 * c * c;
}

static double omega_ai_func_prcorr (CubicEosModel* model, size_t i, double t)
{
  double sqr_tri = sqrt(t/model->tc[i]);
  double fw = model->poly_w[i];
  double oma0 = model->omega_a0[i];
  double c = 1.0 + fw * (1.0 - sqr_tri);
  return oma0 * c * c;
}

static double omega_ai_func_srk (CubicEosModel* model, size_t i, double t)
{
  double sqr_tri = sqrt(t/model->tc[i]);
  double fw = model->poly_w[i];
  double oma0 = model->omega_a0[i];
  double c = 1.0 + fw * (1.0 - sqr_tri);
  return oma0 * c * c;
}

static double omega_ai_func_rk (CubicEosModel* model, size_t i, double t)
{
  double sqr_tri = sqrt(t/model->tc[i]);
  double oma0 = model->omega_a0[i];
  return oma0 / sqr_tri;
}

static void init_cubic_eos_model_pr (CubicEosModel* model)
{
  /* default values for PR */
  model->m1 = 1.0 + sqrt(2.0);
  model->m2 = 1.0 - sqrt(2.0);
  const double omega_a0_default = 0.457235529;
  const double omega_b0_default = 0.07796074;
  for (size_t i=0; i < model->nc; ++i) {
    if (! model->omega_a0[i]) {
      model->omega_a0[i] = omega_a0_default;
    }
    if (! model->omega_b0[i]) {
      model->omega_b0[i] = omega_b0_default;
    }
  }

  /* allocate and calculate poly_w */
  model->poly_w = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "polynomial of eccentric factors");
  for (size_t i=0; i < model->nc; ++i) {
    model->poly_w[i] = compute_poly_w_pr (model->w[i]);
  }
  model->omega_ai_func = omega_ai_func_pr;
}

static void init_cubic_eos_model_prcorr (CubicEosModel* model)
{
  model->m1 = 1.0 + sqrt(2.0);
  model->m2 = 1.0 - sqrt(2.0);
  const double omega_a0_default = 0.457235529;
  const double omega_b0_default = 0.07796074;
  for (size_t i=0; i < model->nc; ++i) {
    if (! model->omega_a0[i]) {
      model->omega_a0[i] = omega_a0_default;
    }
    if (! model->omega_b0[i]) {
      model->omega_b0[i] = omega_b0_default;
    }
  }

  /* allocate and calculate poly_w */
  model->poly_w = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "polynomial of eccentric factors");
  for (size_t i=0; i < model->nc; ++i) {
    model->poly_w[i] = compute_poly_w_prcorr (model->w[i]);
  }
  model->omega_ai_func = omega_ai_func_prcorr;
}

static void init_cubic_eos_model_rk (CubicEosModel* model)
{
  model->m1 = 0.0;
  model->m2 = 1.0;
  const double omega_a0_default = 0.4274802;
  const double omega_b0_default = 0.08664035;
  for (size_t i=0; i < model->nc; ++i) {
    if (! model->omega_a0[i]) {
      model->omega_a0[i] = omega_a0_default;
    }
    if (! model->omega_b0[i]) {
      model->omega_b0[i] = omega_b0_default;
    }
  }

  /* no need of poly_w for RK */
  model->omega_ai_func = omega_ai_func_rk;
}

static void init_cubic_eos_model_srk (CubicEosModel* model)
{
  model->m1 = 0.0;
  model->m2 = 1.0;
  const double omega_a0_default = 0.4274802;
  const double omega_b0_default = 0.08664035;
  for (size_t i=0; i < model->nc; ++i) {
    if (! model->omega_a0[i]) {
      model->omega_a0[i] = omega_a0_default;
    }
    if (! model->omega_b0[i]) {
      model->omega_b0[i] = omega_b0_default;
    }
  }

  /* allocate and calculate poly_w */
  model->poly_w = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "polynomial of eccentric factors");
  for (size_t i=0; i < model->nc; ++i) {
    model->poly_w[i] = compute_poly_w_srk (model->w[i]);
  }
  model->omega_ai_func = omega_ai_func_srk;
}

/* Wrapper for LAPACK function DSPEV */
int dspev (char JOBZ, char UPLO,
           int N, double *A,
           double *W, double *Z,
           int LDZ, double *WORK)
{
  extern void dspev_ (char* JOBZp, char* UPLOp,
                      int* Np, double* A,
                      double* W, double* Z,
                      int* LDZp, double* WORK, int* INFOp);
  int INFO;
  dspev_ (&JOBZ, &UPLO, &N, A, W, Z, &LDZ, WORK, &INFO);
  return INFO;
}

static void update_t(CubicEosModel* model, double t)
{
  for (size_t i=0; i < model->nc; ++i) {
    double tri = t / model->tc[i];
    model->bi_t[i] = model->omega_b0[i] / tri;
    model->sq_ai_t[i] = sqrt(model->omega_ai_func(model, i, t)) / tri;
  }
}

static void update_p(CubicEosModel* model, double p)
{
  for (size_t i=0; i < model->nc; ++i) {
    double pri = p / model->pc[i];
    model->bi[i] = model->bi_t[i] * pri;
    model->sq_ai[i] = model->sq_ai_t[i] * sqrt(pri);
    for (size_t j=0; j < model->rank; ++j) {
      model->q_mat[i*model->nc + j] = model->eigv[j*model->nc + i] * model->sq_ai[i];
    }
    model->q_mat[i*model->nc + model->rank] = model->bi[i];
  }
}

static void update_tp (CubicEosModel* model, double t, double p)
{
  if (t != model->t) {
    model->t = t;
    model->p = p;
    update_t(model, t);
    update_p(model, p);
  } else {
    if (p != model->p) {
      model->p = p;
      update_p(model, p);
    }
  }
}

static void init_cubic_eos_model (CubicEosModel* model)
{
  /* solve eigensystem of D=1-bic only once */
  size_t c = model->nc;
  size_t nd = (size_t) c * (c + 1) / 2;
  double* D = (double*) scm_malloc (nd * sizeof(double));
  /* matrix D is stored as a BLAS packed symmetric matrix using lower triangle and column major */
  for (size_t i=0; i < c; ++i) {
    for (size_t j=0; j <= i; ++j) {
      size_t p = (size_t) (i + (2 * c - j - 1) * j / 2);
      if (j == i) {
        D[p] = 1.0;
      } else {
        size_t q = (size_t) ((i - 1) * i / 2 + j);
        D[p] = 1.0 - model->bic[q];
      }
    }
  }
#ifdef DEBUG
  for (size_t i=0; i < nd; ++i) {
    printf("%6.2f", D[i]);
  }
  putchar('\n');
#endif
  double* W = (double*) scm_malloc (c * sizeof(double));
  double* Z = (double*) scm_malloc (c * c * sizeof(double));
  double* WORK = (double*) scm_malloc (3*c * sizeof(double));
  int info = dspev ('V', 'L', c, D, W, Z, c, WORK);
  if (info != 0) {
    scm_throw (scm_from_utf8_symbol ("error-internal"),
               scm_list_2 (scm_from_utf8_string ("Fail to solve eigen system (dspev return error ~a)"), scm_from_int(info)));
  }
  free (D);
  free (WORK);
#ifdef DEBUG
  for (size_t i=0; i < c; ++i) {
    printf("%15.8f:", W[i]);
    for (size_t j=0; j < c; ++j) {
      printf("%15.8f", Z[i*c+j]);
    }
    putchar('\n');
  }
  putchar('\n');
#endif
  /* Determine which eigenvalues we should keep */
  double norm = cblas_dnrm2 (c, W, 1);
  model->rank = 0;
  size_t* order = (size_t*) scm_malloc (c * sizeof(size_t));
  const double cutoff = scm_to_double(scm_variable_ref(eos_cubic_rv_cutoff));
  for (size_t i=0; i < c; i++) {
    if ( fabs(W[i]/norm) > cutoff ) {
      order[model->rank] = i;
      model->rank++;
    }
  }
  model->lambda = (double*)
    scm_gc_malloc_pointerless (model->rank * sizeof(double), "lambda");
  model->eigv = (double*)
    scm_gc_malloc_pointerless (model->rank * model->nc * sizeof(double), "eigv");
  for (size_t i=0; i < model->rank; ++i) {
    model->lambda[i] = W[order[i]];
    for (size_t j=0; j < c; ++j) {
      model->eigv[i*c+j] = Z[order[i]*c+j];
    }
  }
  free (W);
  free (Z);
#ifdef DEBUG
  printf("RANK=%d\n", model->rank);
  for (size_t i=0; i < model->rank; ++i) {
    printf("%15.8f:", model->lambda[i]);
    for (size_t j=0; j < c; ++j) {
      printf("%15.8f", model->eigv[i*c+j]);
    }
    putchar('\n');
  }
  putchar('\n');
#endif

  /* allocate memory for working zone */
  model->sq_ai_t = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "sq_ai_t");
  model->sq_ai = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "sq_ai");
  model->bi_t = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "bi_t");
  model->bi = (double*)
    scm_gc_malloc_pointerless (model->nc * sizeof(double), "bi");
  model->q_mat = (double*)
    scm_gc_malloc_pointerless (model->nc * model->rank+1 * sizeof(double), "q_mat");

  /* initialize temperature and pressure to standard condition 60F and 14.7psia */
  model->t = 0;            /* rankine */
  model->p = 0;              /* psia */
  update_tp (model, 519.67, 14.7);
#ifdef DEBUG
  for (size_t i=0; i < model->rank+1; ++i) {
    for (size_t j=0; j < c; ++j) {
      printf("%15.8f", model->q_mat[i*model->rank+1+j]);
    }
    putchar('\n');
  }
  putchar('\n');
#endif

}

SCM_DEFINE (eos_cubic_make, "eos-cubic-make", 2, 0, 1,
            (SCM type, SCM nc, SCM rest),
            "Make a cubic EOS model.")
{
  /* parse the keyword arguments */
  SCM k_cnames = scm_from_utf8_keyword ("cnames");
  SCM k_mw = scm_from_utf8_keyword ("mw");
  SCM k_pc = scm_from_utf8_keyword ("pc");
  SCM k_tc = scm_from_utf8_keyword ("tc");
  SCM k_w = scm_from_utf8_keyword ("w");
  SCM k_bic = scm_from_utf8_keyword ("bic");
  SCM k_omega_a0 = scm_from_utf8_keyword ("omega-a");
  SCM k_omega_b0 = scm_from_utf8_keyword ("omega-b");
  SCM k_vshift = scm_from_utf8_keyword ("vshift");

  SCM cnames = SCM_UNDEFINED;
  SCM mw = SCM_UNDEFINED;
  SCM pc = SCM_UNDEFINED;
  SCM tc = SCM_UNDEFINED;
  SCM w = SCM_UNDEFINED;
  SCM bic = SCM_UNDEFINED;
  SCM omega_a0 = SCM_UNDEFINED;
  SCM omega_b0 = SCM_UNDEFINED;
  SCM vshift = SCM_UNDEFINED;

  scm_c_bind_keyword_arguments ("make-cubic-eos-model", rest, 0,
                                k_cnames, &cnames,
                                k_mw, &mw,
                                k_pc, &pc,
                                k_tc, &tc,
                                k_w, &w,
                                k_bic, &bic,
                                k_omega_a0, &omega_a0,
                                k_omega_b0, &omega_b0,
                                k_vshift, &vshift,
                                SCM_UNDEFINED);
  /* test required arguments */
  if (SCM_UNBNDP (cnames) ||
      SCM_UNBNDP (mw) ||
      SCM_UNBNDP (pc) ||
      SCM_UNBNDP (tc) ||
      SCM_UNBNDP (w)) {
    scm_throw (scm_from_utf8_symbol ("error-missing-argument"),
               scm_list_1 (scm_from_utf8_string ("Function make-cubic-eos-model must be called \
with at least the following keyword arguments: #:cnames #:mw #:pc #:tc #:w")));
  }

  struct CubicEosModel* model;
  model = (struct CubicEosModel*)
    scm_gc_malloc (sizeof(struct CubicEosModel), "cubic eos model");
  model->type = type;
  size_t c = scm_to_size_t(nc);
  model->nc = c;
  model->cnames = cnames;

  /* allocate data arrays */
  model->mw = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "molecular weights");
  model->pc = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "critical pressures");
  model->tc = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "critical temperatures");
  model->w = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "accentric factors");
  size_t nbic = (size_t) (c * (c - 1) / 2);
  model->bic = (double*)
    scm_gc_malloc_pointerless (nbic * sizeof(double), "binary interaction coefficients");
  model->omega_a0 = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "OmegaA0s");
  model->omega_b0 = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "OmegaB0s");
  model->vshift = (double*)
    scm_gc_malloc_pointerless (c * sizeof(double), "Volume shifts");

  /* initialize data arrays */
  for (size_t i=0; i < c; ++i) {
    model->mw[i] = scm_to_double (scm_c_vector_ref (mw, i));
    model->pc[i] = scm_to_double (scm_c_vector_ref (pc, i));
    /* convert temperature to rankine */
    model->tc[i] = scm_to_double (scm_c_vector_ref (tc, i)) + 459.67;
    model->w[i] = scm_to_double (scm_c_vector_ref (w, i));
    if (SCM_UNBNDP (omega_a0)) {
      model->omega_a0[i] = 0;
    } else {
      model->omega_a0[i] = scm_to_double (scm_c_vector_ref (omega_a0, i));
    }
    if (SCM_UNBNDP (omega_b0)) {
      model->omega_b0[i] = 0;
    } else {
      model->omega_b0[i] = scm_to_double (scm_c_vector_ref (omega_b0, i));
    }
    if (SCM_UNBNDP (vshift)) {
      model->vshift[i] = 0;
    } else {
      model->vshift[i] = scm_to_double (scm_c_vector_ref (vshift, i));
    }
  }
  for (size_t i=0; i < nbic; ++i) {
    if (SCM_UNBNDP (bic)) {
      model->bic[i] = 0;
    } else {
      model->bic[i] = scm_to_double (scm_c_vector_ref (bic, i));
    }
  }

  /* initialize model according to type */
  if (scm_symbol_is_equal (type, scm_from_utf8_symbol ("PR"))) {
    init_cubic_eos_model_pr (model);
  }
  else if (scm_symbol_is_equal (type, scm_from_utf8_symbol ("PRCORR"))) {
    init_cubic_eos_model_prcorr (model);
  }
  else if (scm_symbol_is_equal (type, scm_from_utf8_symbol ("RK"))) {
    init_cubic_eos_model_rk (model);
  }
  else if (scm_symbol_is_equal (type, scm_from_utf8_symbol ("SRK"))) {
    init_cubic_eos_model_srk (model);
  }
  else {
    scm_throw (scm_from_utf8_symbol ("error-not-implemented"),
               scm_list_2 (scm_from_utf8_string ("User supplied EOS type ~a has to \
be one from PR/PRCORR/RK/SRK."), type));
  }

  /* initialize model for calculation */
  init_cubic_eos_model (model);

  /* return the wraped scm object */
  return scm_make_foreign_object_1 (eos_cubic_model, model);
}

SCM_DEFINE (eos_cubic_print, "eos-cubic-print", 2, 0, 0,
            (SCM port, SCM model),
            "Print cubic EOS model.")
{
  CubicEosModel* m;
  scm_assert_foreign_object_type (eos_cubic_model, model);
  m = scm_foreign_object_ref (model, 0);
  scm_simple_format (port,
                     scm_from_utf8_string ("~a ~a\n"),
                     scm_list_2(m->type, scm_from_size_t(m->nc)));
  for (size_t i=0; i < m->nc; ++i) {
    scm_simple_format (port,
                       scm_from_utf8_string ("~a ~a ~a ~a ~a ~a ~a ~a\n"),
                       scm_list_n (scm_vector_ref(m->cnames, scm_from_size_t(i)),
                                   scm_from_double(m->mw[i]),
                                   scm_from_double(m->pc[i]),
                                   scm_from_double(m->tc[i]),
                                   scm_from_double(m->w[i]),
                                   scm_from_double(m->omega_a0[i]),
                                   scm_from_double(m->omega_b0[i]),
                                   scm_from_double(m->vshift[i]),
                                   SCM_UNDEFINED));
  }
  size_t nbic = (size_t) (m->nc * (m->nc - 1) / 2);
  for (size_t i=0; i < nbic; ++i) {
    scm_simple_format (port,
                       scm_from_utf8_string ("~a "),
                       scm_list_1 (scm_from_double(m->bic[i])));
  }
  scm_simple_format (port, scm_from_utf8_string ("\n"),
                     scm_list_n (SCM_UNDEFINED));
  return SCM_UNDEFINED;
}


void eos_cubic_init (void)
{
  SCM name, slots;
  scm_t_struct_finalize finalizer;

  name = scm_from_utf8_symbol ("eos-cubic-model");
  slots = scm_list_1 (scm_from_utf8_symbol ("data"));
  finalizer = NULL;

  eos_cubic_model =
    scm_make_foreign_object_type (name, slots, finalizer);

#ifndef SCM_MAGIC_SNARFER
#include "eos.x"
#endif
}

