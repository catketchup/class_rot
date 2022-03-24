/** @file rotation.c Documented rotation module
 *
 * Hongbo Cai and Yilun Guan, 09/25/2021
 *
 * This module computes the rotated temperature and polarization
 * anisotropy power spectra \f$ C_l^{X},\f$'s given the
 * unrotated temperature, polarization power spectra, isotropic rotation angle and amplitude of rotation power spectrum.
 *
 * Follows the full-sky method introduced in https://arxiv.org/abs/1303.1881.
 *
 * The following functions can be called from other modules:
 *
 * -# rotation_init() at the beginning (but after harmonic_init())
 * -# rotation_cl_at_l() at any time for computing Cl_rotated at any l
 * -# rotation_free() at the end
 */

#include <stdio.h>
#include "rotation.h"
#include "lensing.h"
#include <time.h>


/**
 * Define rotation power spectrum.
 *
 * @param pro       Input: pointer to rotation structure
 * @param l         Input: multipole number
 * @return the rotation power spectrum at given l
 */

double rotation_cl_aa_at_l(struct rotation * pro, int l) {
  return 2.*_PI_*pro->A_cb/(l*(l+1));
}

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR
 *
 * This routine evaluates all the rotated \f$ C_l\f$'s at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be called from whatever module at whatever time,
 * provided that rotation_init() has been called before, and
 * rotation_free() has not been called yet.
 *
 * @param pro        Input: pointer to rotation structure
 * @param l          Input: multipole number
 * @param cl_rotated  Output: rotated \f$ C_l\f$'s for all types (TT, TE, EE, etc..)
 * @return the error status
 */

int rotation_cl_at_l(
  struct rotation * pro,
  int l,
  double * cl_rotated    /* array with argument cl_rotated[index_ct] (must be already allocated) */
  ) {

  int last_index;
  int index_lt;

  class_test(l > pro->l_rotated_max,
             pro->error_message,
             "you asked for rotated Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",l,pro->l_rotated_max);

  class_call(array_interpolate_spline(pro->l,
                                      pro->l_size,
                                      pro->cl_rot,
                                      pro->ddcl_rot,
                                      pro->lt_size,
                                      l,
                                      &last_index,
                                      cl_rotated,
                                      pro->lt_size,
                                      pro->error_message),
             pro->error_message,
             pro->error_message);

  /* set to zero for the types such that l<l_max */
  for (index_lt=0; index_lt<pro->lt_size; index_lt++)
    if ((int)l > pro->l_max_lt[index_lt])
      cl_rotated[index_lt]=0.;

  return _SUCCESS_;
}

/**
 * This routine initializes the rotation structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppt Input: pointer to perturbation structure (just in case, not used in current version...)
 * @param phr Input: pointer to harmonic structure
 * @param pfo Input: pointer to fourier structure
 * @param pro Output: pointer to initialized rotation structure
 * @return the error status
 */

int rotation_init(
  struct precision * ppr,
  struct perturbations * ppt,
  struct harmonic * phr,
  struct lensing *ple,
  struct fourier * pfo,
  struct rotation * pro
  ) {

  /** Summary: */
  /** - Define local variables */

  double * mu; /* mu[index_mu]: discretized values of mu
                  between -1 and 1, roots of Legendre polynomial */
  double * w8; /* Corresponding Gauss-Legendre quadrature weights */
  double theta,delta_theta;

  double ** d00; /* dmn[index_mu][index_l] */
  double ** d02 = NULL;
  double ** d22 = NULL;
  double ** dm22 = NULL;
  double * buf_dxx; /* buffer */

  double * Ca; /* C_a[index_mu] */
  /* double * W_Ea; /\* W_Ea[index_mu] *\/ */

  double * ksip = NULL; /* ksip[index_mu] */
  double * ksim = NULL; /* ksim[index_mu] */
  double * ksiX = NULL; /* ksiX[index_mu] */

  double * ksip_ptb = NULL; /* ksip[index_mu] */
  double * ksim_ptb = NULL; /* ksim[index_mu] */

  double fac;

  int num_mu,index_mu,icount;
  int l;
  double ll;
  double * cl_unrotated; /* cl_unrotated[index_ct] */
  double * cl_tt = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_te = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_tb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_ee = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_bb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_eb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
  double * cl_aa; /* array to store rotation power spectrum cl_aa */

  double res, resX, rot;
  double resp, resm;

  double resp_ptb, resm_ptb;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;

  /** - check that we really want to compute at least one spectrum */

  if (pro->has_rotated_cls == _FALSE_) {
    if (pro->rotation_verbose > 0)
      printf("No rotation requested. Rotation module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (pro->rotation_verbose > 0) {
      printf("Computing rotated spectra ");
      if (ppr->accurate_rotation==_TRUE_)
        printf("(accurate mode)\n");
      else
        printf("(fast mode)\n");
    }
  }

  /** - initialize indices and allocate some of the arrays in the
      rotation structure */

  class_call(rotation_indices(ppr,phr,pro),
             pro->error_message,
             pro->error_message);

  /** - put all precision variables hare; will be stored later in precision structure */
  /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/

  if (ppr->accurate_rotation == _TRUE_) {
    num_mu=(pro->l_unrotated_max+ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
    num_mu += num_mu%2; /* Force it to be even */
  } else {
    /* Integrate correlation function difference on [0,pi/16] */
    num_mu = (pro->l_unrotated_max * 2 )/16;
  }

  /** - allocate array of \f$ \mu \f$ values, as well as quadrature weights */
  class_alloc(mu,
              num_mu*sizeof(double),
              pro->error_message);
  /* Reserve last element of mu for mu=1, needed for Ca0? */
  mu[num_mu-1] = 1.0;

  class_alloc(w8,
              (num_mu-1)*sizeof(double),
              pro->error_message);

  if (ppr->accurate_rotation == _TRUE_) {

    class_call(quadrature_gauss_legendre(mu,
                                         w8,
                                         num_mu-1,
                                         ppr->tol_gauss_legendre,
                                         pro->error_message),
               pro->error_message,
               pro->error_message);

  } else { /* Crude integration on [0,pi/16]: Riemann sum on theta */
    delta_theta = _PI_/16. / (double)(num_mu-1);
    for (index_mu=0;index_mu<num_mu-1;index_mu++) {
      theta = (index_mu+1)*delta_theta;
      mu[index_mu] = cos(theta);
      w8[index_mu] = sin(theta)*delta_theta; /* We integrate on mu */
    }
  }

  /** - Compute \f$ d^l_{mm'} (\mu) \f$*/

  icount = 0;
  class_alloc(d00,
              num_mu*sizeof(double*),
              pro->error_message);
  icount += num_mu*(pro->l_unrotated_max+1);

  if(pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
    class_alloc(d22,
                num_mu*sizeof(double*),
                pro->error_message);

    class_alloc(dm22,
                num_mu*sizeof(double*),
                pro->error_message);
    icount += 2*num_mu*(pro->l_unrotated_max+1);
  }
  if(pro->has_eb==_TRUE_ && pro->has_ee==_FALSE_ && pro->has_bb == _FALSE_) {
    class_alloc(dm22,
                num_mu*sizeof(double*),
                pro->error_message);
  }

  /** - Allocate main contiguous buffer **/
  class_alloc(buf_dxx,
              icount * sizeof(double),
              pro->error_message);

  icount = 0;
  for (index_mu=0; index_mu<num_mu; index_mu++) {
    d00[index_mu] = &(buf_dxx[icount+index_mu            * (pro->l_unrotated_max+1)]);
  }
  icount +=num_mu*(pro->l_unrotated_max+1);

  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
    for (index_mu=0; index_mu<num_mu; index_mu++){
      d22[index_mu] = &(buf_dxx[icount+index_mu            * (pro->l_unrotated_max+1)]);
      dm22[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)  * (pro->l_unrotated_max+1)]);
    }
    icount +=2*num_mu*(pro->l_unrotated_max+1);
  }

  if (pro->has_eb==_TRUE_ && pro->has_ee==_FALSE_ && pro->has_bb == _FALSE_) {
  	for (index_mu=0; index_mu<num_mu; index_mu++){
      dm22[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)  * (pro->l_unrotated_max+1)]);
  	}
  	icount +=num_mu*(pro->l_unrotated_max+1);
  }


  class_call(lensing_d00(mu,num_mu,pro->l_unrotated_max,d00),
             pro->error_message,
             pro->error_message);

  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
    class_call(lensing_d22(mu,num_mu,pro->l_unrotated_max,d22),
               pro->error_message,
               pro->error_message);
    class_call(lensing_d2m2(mu,num_mu,pro->l_unrotated_max,dm22),
               pro->error_message,
               pro->error_message);
  }

  if (pro->has_eb==_TRUE_ && pro->has_bb==_FALSE_ && pro->has_ee==_FALSE_) {
    class_call(lensing_d2m2(mu,num_mu,pro->l_unrotated_max,dm22),
               pro->error_message,
               pro->error_message);
  }

  /** - compute \f$ Ca(\mu)\f$ */
  class_alloc(Ca,
              num_mu*sizeof(double),
              pro->error_message);

  class_alloc(cl_unrotated,
              phr->ct_size*sizeof(double),
              pro->error_message);

  /** - Locally store unrotated temperature \f$ cl\f$ and potential \f$ cl_{aa}\f$ spectra **/
  class_alloc(cl_tt,
              (pro->l_unrotated_max+1)*sizeof(double),
              pro->error_message);
  if (pro->has_te==_TRUE_) {
    class_alloc(cl_te,
                (pro->l_unrotated_max+1)*sizeof(double),
                pro->error_message);
  }
  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
    class_alloc(cl_ee,
                (pro->l_unrotated_max+1)*sizeof(double),
                pro->error_message);
    class_alloc(cl_bb,
                (pro->l_unrotated_max+1)*sizeof(double),
                pro->error_message);
  }
  class_alloc(cl_aa,
              (pro->l_unrotated_max+1)*sizeof(double),
              pro->error_message);

  class_alloc(cl_md_ic,
              phr->md_size*sizeof(double *),
              pro->error_message);

  class_alloc(cl_md,
              phr->md_size*sizeof(double *),
              pro->error_message);

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)

      class_alloc(cl_md[index_md],
                  phr->ct_size*sizeof(double),
                  pro->error_message);

    if (phr->ic_size[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                  pro->error_message);
  }

  /* use input cl_aa or generate scale-invariant cl_aa */
  if (pro->claa_from_file==_TRUE_) {
    FILE *f;
    f = fopen(pro->input_claa, "r");
    for (l=0; l<=pro->l_unrotated_max; l++) {
      fscanf(f, "%lf", &cl_aa[l]);
    }
    fclose(f);
  }
  else{
    for (l=1; l<=pro->l_unrotated_max; l++) {
      cl_aa[l] = rotation_cl_aa_at_l(pro, l);
  }
  }

  /* compute rotated lensed Cl, if "lensing" is given */
  if (ppt->has_cl_cmb_lensing_potential == _TRUE_){
    /* note here we l_lensed_max = l_rotated_max, we have to use l<=pro->l_rotated_max, so the last num_mu_minus_lmax numbers are not accurate */
    for (l=2; l<=pro->l_rotated_max; l++) {
      class_call(lensing_cl_at_l(ple,l,cl_unrotated),
               ple->error_message,
               pro->error_message);
    }

    cl_tt[l] = cl_unrotated[pro->index_lt_tt];
    if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
      cl_ee[l] = cl_unrotated[pro->index_lt_ee];
      cl_bb[l] = cl_unrotated[pro->index_lt_bb];
    }
    if (pro->has_te==_TRUE_) {
      cl_te[l] = cl_unrotated[pro->index_lt_te];
    }
  }

  /* else, compute rotated raw cl */
  else {
    for (l=2; l<=pro->l_unrotated_max; l++) {
      class_call(harmonic_cl_at_l(phr,l,cl_unrotated,cl_md,cl_md_ic),
                 phr->error_message,
                 pro->error_message);
    }

    cl_tt[l] = cl_unrotated[pro->index_lt_tt];
    if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
      cl_ee[l] = cl_unrotated[pro->index_lt_ee];
      cl_bb[l] = cl_unrotated[pro->index_lt_bb];
    }
    if (pro->has_te==_TRUE_) {
      cl_te[l] = cl_unrotated[pro->index_lt_te];
    }
  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)
      free(cl_md[index_md]);

    if (phr->ic_size[index_md] > 1)
      free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /** - Compute Ca\f$(\mu)\f$ **/

#pragma omp parallel for                        \
  private (index_mu,l)                          \
  schedule (static)
  for (index_mu=0; index_mu<num_mu; index_mu++) {
    Ca[index_mu] = 0;

    for (l=1; l<=pro->l_unrotated_max; l++) {
      Ca[index_mu] += (2.*l+1.)*cl_aa[l]*d00[index_mu][l];
    }

    Ca[index_mu] /= 4.*_PI_;
  }

  /** - compute ksi+, ksi-, ksiX for rotation*/

  /** - --> ksip, ksim for EE, BB **/
  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {

    class_calloc(ksip,
                 (num_mu-1),
                 sizeof(double),
                 pro->error_message);

    class_calloc(ksim,
                 (num_mu-1),
                 sizeof(double),
                 pro->error_message);

    if (pro->perturb_rotation==_TRUE_){

      class_calloc(ksip_ptb,
                   (num_mu-1),
                   sizeof(double),
                   pro->error_message);
      class_calloc(ksim_ptb,
                   (num_mu-1),
                   sizeof(double),
                   pro->error_message);
    }

  }

  /** - --> ksiX is for EB **/
  if (pro->has_eb==_TRUE_) {
    class_calloc(ksiX,
                 (num_mu-1),
                 sizeof(double),
                 pro->error_message);
  }

#pragma omp parallel for                            \
  private (index_mu, l, ll, resp, resm, resX, fac)	\
  schedule (static)

  for (index_mu=0;index_mu<num_mu-1;index_mu++) {

    for (l=2;l<=pro->l_unrotated_max;l++) {

      ll = (double)l;

      fac = 2*ll+1;

      if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {

        resp = fac*d22[index_mu][l]*(cl_ee[l]+cl_bb[l]);
        resm = fac*dm22[index_mu][l]*(cl_ee[l]-cl_bb[l]);

        ksip[index_mu] += resp;
        ksim[index_mu] += resm;

        if (pro->perturb_rotation==_TRUE_){
          resp_ptb = fac*d22[index_mu][l]*cl_ee[l];
          resm_ptb = fac*dm22[index_mu][l]*cl_ee[l];

          ksip_ptb[index_mu] += resp_ptb;
          ksim_ptb[index_mu] += resm_ptb;
        }
      }

      if (pro->has_eb==_TRUE_) {
        resX = fac*dm22[index_mu][l]*(cl_ee[l]-cl_bb[l]);

        ksiX[index_mu] += resX;
      }
    }

    ksip[index_mu] *= exp(4*Ca[index_mu]);
    ksim[index_mu] *= exp(-4*Ca[index_mu]);

    if (pro->perturb_rotation==_TRUE_) {
      ksip_ptb[index_mu] *= Ca[index_mu];
      ksim_ptb[index_mu] *= Ca[index_mu];
    }
    if (pro->has_eb == _TRUE_) {
      ksiX[index_mu] *= exp(-4*Ca[index_mu]);
    }
  }


  /** - compute rotated \f$ C_l\f$'s by multiplation or integration */
  if (pro->has_tt==_TRUE_) {
    class_call(rotation_rotated_cl_tt(cl_tt,pro),
               pro->error_message,
               pro->error_message);

  }

  if (pro->has_te==_TRUE_) {
    class_call(rotation_rotated_cl_te(cl_te,pro->alpha,Ca[0],pro),
               pro->error_message,
               pro->error_message);
  }


  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
    class_call(rotation_rotated_cl_ee_bb(ksip,ksim,d22,dm22,w8,pro->alpha,Ca[0],num_mu-1,pro),
               pro->error_message,
               pro->error_message);
    if (pro->perturb_rotation==_TRUE_) {
      class_call(rotation_rotated_cl_bb_perturb(ksip_ptb,ksim_ptb,d22,dm22,w8,num_mu-1,pro),
                 pro->error_message,
                 pro->error_message);
    }

  }

  if (pro->has_tb==_TRUE_) {
    pro->l_max_lt[pro->index_lt_tb] = ppt->l_scalar_max;
    class_call(rotation_rotated_cl_tb(cl_te,pro->alpha,Ca[0],pro),
               pro->error_message,
               pro->error_message);
  }

  if (pro->has_eb==_TRUE_) {
    pro->l_max_lt[pro->index_lt_eb] = ppt->l_scalar_max;
    class_call(rotation_rotated_cl_eb(ksiX,dm22,w8,pro->alpha,Ca[0],num_mu-1,pro),
               pro->error_message,
               pro->error_message);
  }

  /** - spline computed \f$ C_l\f$'s in view of interpolation */

  class_call(array_spline_table_lines(pro->l,
                                      pro->l_size,
                                      pro->cl_rot,
                                      pro->lt_size,
                                      pro->ddcl_rot,
                                      _SPLINE_EST_DERIV_,
                                      pro->error_message),
             pro->error_message,
             pro->error_message);

  /** - Free lots of stuff **/
  free(buf_dxx);
  free(d00);
  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_ || pro->has_eb==_TRUE_) {
    free(d22);
    free(dm22);
  }
  free(Ca);
  free(mu);
  free(w8);

  free(cl_unrotated);
  free(cl_tt);
  if (pro->has_te==_TRUE_ || pro->has_tb==_TRUE_)
    free(cl_te);
  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_ || pro->has_eb==_TRUE_) {
    free(cl_ee);
    free(cl_bb);
  }
  free(cl_aa);
  /** - Exit **/

  return _SUCCESS_;

}

/**
 * This routine frees all the memory space allocated by rotation_init().
 *
 * To be called at the end of each run, only when no further calls to
 * rotation_cl_at_l() are needed.
 *
 * @param pro Input: pointer to rotation structure (which fields must be freed)
 * @return the error status
 */

int rotation_free(
  struct rotation * pro
  ) {

  if (pro->has_rotated_cls == _TRUE_) {

    free(pro->l);
    free(pro->cl_rot);
    free(pro->ddcl_rot);
    free(pro->l_max_lt);

  }

  return _SUCCESS_;

}

int rotation_indices(
  struct precision * ppr,
  struct harmonic * phr,
  struct rotation * pro
  ){

  int index_l;

  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */

  int index_md;
  int index_ct;

  /* indices of all Cl types (rotated and unrotated) */

  if (phr->has_tt == _TRUE_) {
    pro->has_tt = _TRUE_;
    pro->index_lt_tt=phr->index_ct_tt;
  }
  else {
    pro->has_tt = _FALSE_;
  }

  if (phr->has_ee == _TRUE_) {
    pro->has_ee = _TRUE_;
    pro->index_lt_ee=phr->index_ct_ee;
  }
  else {
    pro->has_ee = _FALSE_;
  }

  if (phr->has_te == _TRUE_) {
    pro->has_te = _TRUE_;
    pro->index_lt_te=phr->index_ct_te;
  }
  else {
    pro->has_te = _FALSE_;
  }

  if (phr->has_bb == _TRUE_) {
    pro->has_bb = _TRUE_;
    pro->index_lt_bb=phr->index_ct_bb;
  }
  else {
    pro->has_bb = _FALSE_;
  }

  if (phr->has_tb == _TRUE_) {
    pro->has_tb = _TRUE_;
    pro->index_lt_tb=phr->index_ct_tb;
  }

  if (phr->has_eb == _TRUE_) {
    pro->has_eb = _TRUE_;
    pro->index_lt_eb=phr->index_ct_eb;
  }

  pro->lt_size = phr->ct_size;

  /* number of multipoles */

  pro->l_unrotated_max = phr->l_max_tot;

  pro->l_rotated_max = pro->l_unrotated_max - ppr->delta_l_max;

  for (index_l=0; (index_l < phr->l_size_max) && (phr->l[index_l] <= pro->l_rotated_max);index_l++);

  if (index_l < phr->l_size_max) index_l++; /* one more point in order to be able to interpolate till pro->l_rotated_max */

  pro->l_size = index_l+1;

  class_alloc(pro->l, pro->l_size*sizeof(double),pro->error_message);

  for (index_l=0; index_l < pro->l_size; index_l++) {

    pro->l[index_l] = phr->l[index_l];

  }

  /* allocate table where results will be stored */

  class_alloc(pro->cl_rot,
              pro->l_size*pro->lt_size*sizeof(double),
              pro->error_message);

  class_alloc(pro->ddcl_rot,
              pro->l_size*pro->lt_size*sizeof(double),
              pro->error_message);

  /* fill with unrotated cls */

  class_alloc(cl_md_ic,
              phr->md_size*sizeof(double *),
              pro->error_message);

  class_alloc(cl_md,
              phr->md_size*sizeof(double *),
              pro->error_message);

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)

      class_alloc(cl_md[index_md],
                  phr->ct_size*sizeof(double),
                  pro->error_message);

    if (phr->ic_size[index_md] > 1)

      class_alloc(cl_md_ic[index_md],
                  phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                  pro->error_message);
  }

  for (index_l=0; index_l<pro->l_size; index_l++) {

    class_call(harmonic_cl_at_l(phr,pro->l[index_l],&(pro->cl_rot[index_l*pro->lt_size]),cl_md,cl_md_ic),
               phr->error_message,
               pro->error_message);

  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)
      free(cl_md[index_md]);

    if (phr->ic_size[index_md] > 1)
      free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /* we want to output Cl_rotated up to the same l_max as Cl_unroated
     (even if a number delta_l_max of extra values of l have been used
     internally for more accurate results). Notable exception to the
     above rule: ClBB_rotated(scalars) must be outputed at least up to the same l_max as
     ClEE_unrotated(scalars) (since ClBB_unroated is null for scalars)
  */

  class_alloc(pro->l_max_lt,pro->lt_size*sizeof(double),pro->error_message);
  for (index_ct = 0; index_ct < pro->lt_size; index_ct++) {
    pro->l_max_lt[index_ct]=0.;
    for (index_md = 0; index_md < phr->md_size; index_md++) {
      pro->l_max_lt[index_ct]=MAX(pro->l_max_lt[index_ct],phr->l_max_ct[index_md][index_ct]);

      if ((pro->has_bb == _TRUE_) && (pro->has_ee == _TRUE_) && (index_ct == pro->index_lt_bb)) {
        pro->l_max_lt[index_ct]=MAX(pro->l_max_lt[index_ct],phr->l_max_ct[index_md][pro->index_lt_ee]);
      }

    }
  }

  return _SUCCESS_;

}

/**
 * rotation does not change Cl_TT
 * @param cl_tt  Input: unrotated cl_tt
 * @param pro   Input/output: Pointer to the rotation structure
 */
int rotation_rotated_cl_tt(double *cl_tt,
                           struct rotation * pro){
  int index_l;
  for(index_l=0; index_l<pro->l_size; index_l++){
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_tt] = cl_tt[(int)pro->l[index_l]];
  }

  return _SUCCESS_;
}

/**
 * This routine computes the rotation power spectra by Gaussian quadrature
 *
 * @param cl_te  Input: unrotated cl_te
 * @param alpha Input: Isotropic rotation angle
 * @param Ca0   Input:
 * @param nmu   Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro   Input/output: Pointer to the rotation structure
 * @return the error status
 */
int rotation_rotated_cl_te(double *cl_te,
                           double alpha,
                           double Ca0,
                           struct rotation * pro
  ){
  int index_l;
  for(index_l=0; index_l<pro->l_size; index_l++){
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_te] = cl_te[(int)pro->l[index_l]]*cos(2*pro->alpha)*exp(-2*Ca0);
  }

  return _SUCCESS_;
}

/**
 * This routine computes the rotated power spectra by Gaussian quadrature
 *
 * @param cl_te  Input: unrotated cl_te
 * @param alpha Input: Isotropic rotation angle
 * @param Ca0   Input:
 * @param nmu   Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro   Input/output: Pointer to the rotation structure
 * @return the error status
 */
int rotation_rotated_cl_tb(double *cl_te,
                           double alpha,
                           double Ca0,
                           struct rotation * pro
  ){
  int index_l;
  for(index_l=0; index_l<pro->l_size; index_l++){
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_tb] = cl_te[(int)pro->l[index_l]]*sin(2*pro->alpha)*exp(-2*Ca0);
  }

  return _SUCCESS_;
}

/**
 * This routine computes the rotation power spectra by Gaussian quadrature
 *
 * @param ksip  Input: rotated correlation function (ksi+[index_mu])
 * @param ksim  Input: rotated correlation function (ksi-[index_mu])
 * @param d22   Input: Wigner d-function (\f$ d^l_{22}\f$[l][index_mu])
 * @param dm22  Input: Wigner d-function (\f$ d^l_{-22}\f$[l][index_mu])
 * @param w8    Input: Legendre quadrature weights (w8[index_mu])
 * @param alpha Input: Isotropic rotation angle
 * @param Ca0   Input: (C^{\alapha{0}})
 * @param nmu   Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro   Input/output: Pointer to the rotation structure
 * @return the error status
 */
int rotation_rotated_cl_ee_bb(double *ksip,
                              double *ksim,
                              double **d22,
                              double **dm22,
                              double *w8,
                              double alpha,
                              double Ca0,
                              int nmu,
                              struct rotation * pro
  ){

  double clp, clm;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clp,clm)                 \
  schedule (static)

  for(index_l=0; index_l < pro->l_size; index_l++){
    clp=0; clm=0;
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][(int)pro->l[index_l]]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*dm22[imu][(int)pro->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    clp *= exp(-4*Ca0)/2;
    clm *= exp(-4*Ca0)*cos(4*alpha)/2;
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_ee]=(clp+clm)/2;
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_bb]=(clp-clm)/2;
  }

  return _SUCCESS_;
}

/**
 * This routine computes the rotated power spectra by Gaussian quadrature
 *
 * @param ksiX Input: Rotated correlation function (ksiX[index_mu])
 * @param dm22  Input: Wigner d-function (\f$ d^l_{-22}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro  Input/output: Pointer to the rotation structure
 * @return the error status
 */
int rotation_rotated_cl_eb(double *ksiX,
                           double **dm22,
                           double *w8,
                           double alpha,
                           double Ca0,
                           int nmu,
                           struct rotation * pro
  ){

  double clX;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clX)                     \
  schedule (static)

  for(index_l=0; index_l < pro->l_size; index_l++){
    clX=0;
    for (imu=0;imu<nmu;imu++) {
      clX += ksiX[imu]*dm22[imu][(int)pro->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_eb]=clX*exp(-4*Ca0)*sin(4*alpha)/4;
  }

  return _SUCCESS_;
}

/**
 * This routine computes the rotation power spectra by Gaussian quadrature using perturbative method
 *
 * @param ksip_ptb  Input:
 * @param ksim_ptb  Input:
 * @param d22   Input: Wigner d-function (\f$ d^l_{22}\f$[l][index_mu])
 * @param dm22  Input: Wigner d-function (\f$ d^l_{-22}\f$[l][index_mu])
 * @param w8    Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu   Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro   Input/output: Pointer to the rotation structure
 * @return the error status
 */

int rotation_rotated_cl_bb_perturb(double *ksip_ptb,
                                   double *ksim_ptb,
                                   double **d22,
                                   double **dm22,
                                   double *w8,
                                   int nmu,
                                   struct rotation * pro
  ){
  double clbb_ptb;
  int imu;
  int index_l;

  /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
  private (imu,index_l,clbb_ptb)                \
  schedule (static)

  for(index_l=0; index_l < pro->l_size; index_l++){
    clbb_ptb=0;
    for (imu=0;imu<nmu;imu++){
      clbb_ptb += (ksip_ptb[imu]*d22[imu][(int)pro->l[index_l]] + ksim_ptb[imu]*dm22[imu][(int)pro->l[index_l]])*w8[imu];
    }
    pro->cl_rot[index_l*pro->lt_size+pro->index_lt_bb]=clbb_ptb;
  }

  return _SUCCESS_;
}
