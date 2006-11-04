/******* d_congrad5_fn_qop_two_src.c - conjugate gradient for SU3/fermions **/
/* MIMD version 7 */

/* This is the two-source MILC wrapper for the SciDAC Level 3 QOP inverter 
   using the raw interface */
/* 2/2005 D. Renner and C. Jung */
/* 5/2005 C. DeTar two source version eliminates one remapping */
/* 9/2005 C. DeTar converted to C code */
/* 12/2005 C. DeTar upgraded to new Level 3 API */

#include "generic_ks_includes.h"
#include <qop.h>

/* The standard MILC interface for two sources */

int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset milc_src1,     /* source vector (type su3_vector) */
    field_offset milc_src2,
    field_offset milc_sol1,	/* solution vectors */
    field_offset milc_sol2,
    Real mass1,
    Real mass2,
    int niter,		        /* maximal number of CG interations */
    int nrestart,               /* maximal number of CG restarts */
    Real rsqmin,	        /* desired residue squared */
    int milc_parity,		/* parity to be worked on */
    Real  *final_rsq_ptr 	/* final residue squared */
    )
{
  int iterations_used;
  static Real t_mass1;
  static Real t_mass2;
  Real *masses[2];
  int nmass[2], nsrc;
  field_offset milc_srcs[2], milc_sols0[1], milc_sols1[1], *milc_sols[2];

  /* Set up general source and solution pointers for two sources,
     one mass per source so one solution per source */
  nsrc = 2;
  milc_srcs[0] = milc_src1;
  milc_srcs[1] = milc_src2;

  nmass[0] = 1; 
  nmass[1] = 1;
  t_mass1 = mass1;
  t_mass2 = mass2;
  masses[0] = &t_mass1;
  masses[1] = &t_mass2;

  milc_sols0[0] = milc_sol1;
  milc_sols1[0] = milc_sol2;
  milc_sols[0] = milc_sols0;
  milc_sols[1] = milc_sols1;

  iterations_used = ks_congrad_qop_site2site( niter, nrestart, rsqmin, 
					      masses, nmass, milc_srcs,
					      milc_sols, nsrc, final_rsq_ptr,
					      milc_parity );

  total_iters += iterations_used;
  return( iterations_used );
}
