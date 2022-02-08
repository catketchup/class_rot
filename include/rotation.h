/** @file rotation.h Documented includes for harmonic module */

#ifndef __ROTATION__
#define __ROTATION__

#include "harmonic.h"
#include "lensing.h"

/**
 * Structure containing everything about rotated spectra that other modules need to know.
 *
 * Once initialized by rotation_init(), contains a table of all rotated
 * \f$ C_l\f$'s for the all modes (scalar/tensor), all types (TT, TE...),
 * and all pairs of initial conditions (adiabatic, isocurvatures...).
 * FOR THE MOMENT, ASSUME ONLY SCALAR & ADIABATIC
 */

struct rotation {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{

  short has_rotated_cls; /**< do we need to compute rotated \f$ C_l\f$'s at all ? */
  short has_cl_cmb_rotation_spectrum; /**< do we need \f$ C_l \f$'s for CMB rotation power spectrum? */

  //@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{
  int use_lensed; /* use lensed unrotated C_l or unlensed unrotated C_l */

  double alpha;
  double A_cb;

  int has_tt;
  int has_te;
  int has_tb;
  int has_ee;
  int has_bb;
  int has_eb;
  int has_aa;
  int claa_from_file; /**do we use input claa txt file */
  FileName input_claa; /**input claa txt file */
  int has_ea;
  int perturb_rotation; /**< do we use the perturbative method */

  int index_lt_tt; /**< index for type \f$ C_l^{TT} \f$*/
  int index_lt_te; /**< index for type \f$ C_l^{TE} \f$*/
  int index_lt_tb; /**< index for type \f$ C_l^{TB} \f$*/
  int index_lt_ee; /**< index for type \f$ C_l^{EE} \f$*/
  int index_lt_bb; /**< index for type \f$ C_l^{BB} \f$*/
  int index_lt_eb; /**< index for type \f$ C_l^{EB} \f$*/
  int index_lt_aa; /**< index for type \f$ C_l^{aa} \f$*/
  int index_lt_ea; /**< index for type \f$ C_l^{ea} \f$*/

  int lt_size; /**< number of \f$ C_l\f$ types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{

  int l_unrotated_max; /**< last multipole in all calculations (same as in harmonic module)*/

  int l_rotated_max; /**< last multipole at which rotated spectra are computed */

  /* interpolable version: */

  int l_size;       /**< number of l values */

  int * l_max_lt;    /**< last multipole (given as an input) at which
                        we want to output \f$ C_l \f$'s for a given mode and type */

  double * l;       /**< table of multipole values l[index_l] */
  double * cl_rot; /**< table of anisotropy spectra for each
                      multipole and types,
                      cl[index_l * ple->lt_size + index_lt] */

  double * ddcl_rot; /**< second derivatives for interpolation */

  //@}

  /** @name - technical parameters */

  //@{

  short rotation_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};


#ifdef __cplusplus
extern "C" {
#endif

  int rotation_cl_at_l(struct rotation * pro,
                       int l,
                       double * cl_rotated);

  int rotation_init(struct precision * ppr,
                    struct perturbations * ppt,
                    struct harmonic * phr,
                    struct fourier * pfo,
                    struct rotation * pro);

  int rotation_free(struct rotation * pro);

  int rotation_indices(struct precision * ppr,
                       struct harmonic * phr,
                       struct rotation * pro);

  /* int rotation_rotated_cl_tt( */
  /* 	); */
  int rotation_rotated_cl_tt(double *cl_tt,
                             struct rotation * pro);

  int rotation_rotated_cl_te(double *cl_te,
                             double alpha,
                             double Ca0,
                             struct rotation * pro);

  int rotation_rotated_cl_tb(double *cl_te,
                             double alpha,
                             double Ca0,
                             struct rotation * pro);

  int rotation_rotated_cl_ee_bb(double *ksip,
                                double *ksim,
                                double **d22,
                                double **dm22,
                                double *w8,
                                double alpha,
                                double Ca0,
                                int nmu,
                                struct rotation * pro);


  int rotation_rotated_cl_eb(double *ksiX,
                             double **dm22,
                             double *w8,
                             double alpha,
                             double Ca0,
                             int nmu,
                             struct rotation * pro);


  int rotation_rotated_cl_bb_perturb(double *ksip_ptb,
                                     double *ksim_ptb,
                                     double **d22,
                                     double **dm22,
                                     double *w8,
                                     int nmu,
                                     struct rotation * pro);
#ifdef __cplusplus
}
#endif

#endif
