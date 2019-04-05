/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//  $Log: setup.c,v $
*/
/* MIMD version 7 */
#define IF_OK if(status==0)
#define _POSIX_C_SOURCE 200112L // for gethostname

#include "ks_imp_includes.h"
#include "quark_action.h"
#include "lattice_qdp.h"
#include "../include/fermion_links.h"
#define SU3_MAT_OP_NO_STORAGE
#include "../include/su3_mat_op.h"
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

/* Each node has a params structure for passing simulation parameters */
#include "params.h"

/* Forward declarations */
static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers(void);

int setup(void)
{
  printf("Amit MyFFApp/setup.c start of function int setup(void) \n");
  int prompt;
  
  /* print banner, get volume, seed */
  prompt = initial_set();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);

  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  node0_printf("Made 3nn gathers\n"); fflush(stdout);
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  
  node0_printf("Finished setup\n"); fflush(stdout);
  printf("Amit MyFFApp/setup.c end of function int setup(void) \n");
  return( prompt );
}

static int n_naiks;
static double eps_naik[MAX_NAIK];

/* SETUP ROUTINES */
static int initial_set(void)
{
  printf("Amit MyFFApp/setup.c start of function static int initial_set(void)  \n");
  int prompt,status,i,tmporder;
  Real current_naik_epsilon;

  /* On node zero, read lattice size, seed, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("Rational function hybrid Monte Carlo algorithm\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");

    /* Print list of options selected */
    node0_printf("Options selected...\n");
    show_generic_opts();
    show_generic_ks_opts();
    show_generic_ks_md_opts();
#ifdef INT_ALG
    node0_printf("INT_ALG=%s\n",ks_int_alg_opt_chr());
#endif
#if FERM_ACTION == HISQ
    show_su3_mat_opts();
    show_hisq_links_opts();
    show_hisq_force_opts();
#elif FERM_ACTION == HYPISQ
    show_su3_mat_opts();
    show_hypisq_links_opts();
    show_hypisq_force_opts();
#endif

    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   par_buf.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   par_buf.ionode_geometry, 4);
#endif
#endif
    IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );
    /* Number of pseudofermions */
    IF_OK status += get_i(stdin, prompt,"n_pseudo", &par_buf.n_pseudo );
    if(par_buf.n_pseudo > MAX_N_PSEUDO){
      printf("Error:  Too many pseudofermion fields.  Recompile. Current max is %d\n"
	     ,MAX_N_PSEUDO);
      terminate(1);
    }
    /* get name of file containing rational function parameters */
    IF_OK status += get_s(stdin, prompt, "load_rhmc_params", 
			  par_buf.rparamfile);
    /* beta, quark masses */
    IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

    IF_OK status += get_i(stdin, prompt,"n_dyn_masses", &par_buf.n_dyn_masses );
    IF_OK status += get_vf(stdin, prompt, "dyn_mass", par_buf.dyn_mass, par_buf.n_dyn_masses);
    IF_OK status += get_vi(stdin, prompt, "dyn_flavors", par_buf.dyn_flavors, par_buf.n_dyn_masses);

    IF_OK status += get_f(stdin, prompt,"u0", &par_buf.u0 );

    if(status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(mynode()==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )
    normal_exit(0);
  
  nx        = par_buf.nx;
  ny        = par_buf.ny;
  nz        = par_buf.nz;
  nt        = par_buf.nt;
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = par_buf.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = par_buf.ionode_geometry[i];
#endif
#endif
  iseed     = par_buf.iseed;
  n_pseudo  = par_buf.n_pseudo;
  strcpy(rparamfile,par_buf.rparamfile);
  
  this_node = mynode();
  number_of_nodes = numnodes();
  volume=nx*ny*nz*nt;
  total_iters=0;
#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = 0;
#endif
      
#ifdef HYPISQ_SVD_COUNTER
  hypisq_svd_counter = 0;
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = 0;
#endif

#ifdef HYPISQ_FORCE_FILTER_COUNTER
  hypisq_force_filter_counter = 0;
#endif

  /* Load rational function parameters */
  rparam = load_rhmc_params(rparamfile, n_pseudo);  
  if(rparam == NULL)terminate(1);

  /* Determine the maximum rational fcn order */
  max_rat_order = 0;
  for(i = 0; i < n_pseudo; i++){
    if(rparam[i].MD.order > max_rat_order)max_rat_order = rparam[i].MD.order;
    if(rparam[i].GR.order > max_rat_order)max_rat_order = rparam[i].GR.order;
    if(rparam[i].FA.order > max_rat_order)max_rat_order = rparam[i].FA.order;
  }
  node0_printf("Maximum rational func order is %d\n",max_rat_order);

  /* Determine the number of different Naik masses
     and fill in n_orders_naik and n_pseudo_naik        */
  current_naik_epsilon = rparam[0].naik_term_epsilon;
  tmporder = 0;
  n_naiks = 0;
  n_order_naik_total = 0;
  for( i=0; i<n_pseudo; i++ ) {
    if( rparam[i].naik_term_epsilon != current_naik_epsilon ) {
      if( tmporder > 0 ) {
        n_orders_naik[n_naiks] = tmporder;
	eps_naik[n_naiks] = current_naik_epsilon;
        current_naik_epsilon = rparam[i].naik_term_epsilon;
        n_naiks++;
        n_order_naik_total += tmporder;
        tmporder = 0;
      }
    }
    tmporder += rparam[i].MD.order;
    n_pseudo_naik[n_naiks]++;
  }
  if( tmporder > 0 ) {
    n_orders_naik[n_naiks] = tmporder;
    eps_naik[n_naiks] = current_naik_epsilon;
    n_order_naik_total += tmporder;
    n_naiks++;
  }
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  // calculate epsilon corrections for different Naik terms
  if( 0 != eps_naik[0] ) {
    node0_printf("IN HISQ AND HYPISQ ACTIONS FIRST SET OF PSEUDO FERMION FIELDS SHOULD HAVE EPSILON CORRECTION TO NAIK TERM ZERO.\n");
    terminate(1);
  }
#endif
  node0_printf("Naik term correction structure of multi_x:\n");
  node0_printf("n_naiks %d\n",n_naiks);
  for( i=0; i<n_naiks; i++ ) {
    node0_printf("n_pseudo_naik[%d]=%d\n", i, n_pseudo_naik[i]);
    node0_printf("n_orders_naik[%d]=%d\n", i, n_orders_naik[i]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
    node0_printf("eps_naik[%d]=%f\n", i, eps_naik[i]);
#endif
  }
  node0_printf("n_order_naik_total %d\n",n_order_naik_total);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  if( n_naiks+1 > MAX_NAIK ) {
    node0_printf("MAX_NAIK=%d < n_naiks+1=%d\n", MAX_NAIK, n_naiks+1 );
    node0_printf("Increase MAX_NAIK\n");
    terminate(1);
  }
#else /* non HISQ */
  if( n_naiks>1 ) {
    node0_printf("FOR ACTIONS OTHER THAN HISQ AND HYPISQ EPSILON CORRECTION IS NOT USED.\n");
    node0_printf("ONLY ONE SET OF X LINKS IS USED.\n");
    node0_printf("SET ALL naik_mass TO 0 IN RATIONAL FUNCTION FILE.\n");
    terminate(1);
  }
#endif /* HISQ */

  beta = par_buf.beta;
  
  n_dyn_masses = par_buf.n_dyn_masses;
  for(i = 0; i < n_dyn_masses; i++){
    dyn_mass[i] = par_buf.dyn_mass[i];
    dyn_flavors[i] = par_buf.dyn_flavors[i];
  }
  u0 = par_buf.u0;
  printf("Amit MyFFApp/setup.c end of function static int initial_set(void)  \n");
  return(prompt);
}

/* read in parameters and coupling constants	*/
int readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  printf("Amit MyFFApp/setup.c start of function int readin(int prompt) \n");
  int status;
  int i;
  
  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {
    
    printf("\n\n");
    status=0;
    
    /* warms, trajecs */
    IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
    IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );
    
    /* trajectories between propagator measurements */
    IF_OK status += 
      get_i(stdin, prompt,"traj_between_meas", &par_buf.propinterval );
    
    /* microcanonical time step */
    IF_OK status += 
      get_f(stdin, prompt,"microcanonical_time_step", &par_buf.epsilon );
    
    /*microcanonical steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );
    
    /* Data for each pseudofermion */

    for(i = 0; i < par_buf.n_pseudo; i++){
      Real tmp[3]; int itmp[3];

      /* Residuals for multicg solves */
      IF_OK status += get_vf(stdin, prompt,"cgresid_md_fa_gr", tmp, 3 );
      /* rsqmin is r**2 in conjugate gradient */
      IF_OK {
	par_buf.rsqmin_md[i] = tmp[0]*tmp[0];
	par_buf.rsqmin_fa[i] = tmp[1]*tmp[1];
	par_buf.rsqmin_gr[i] = tmp[2]*tmp[2];
      }

      /* Max CG iterations for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "max_multicg_md_fa_gr", itmp, 3);
      IF_OK {
	par_buf.niter_md[i] = itmp[0];
	par_buf.niter_fa[i] = itmp[1];
	par_buf.niter_gr[i] = itmp[2];
      }

      /* Precision for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "cgprec_md_fa_gr", itmp, 3);
      IF_OK {
	par_buf.prec_md[i] = itmp[0];
	par_buf.prec_fa[i] = itmp[1];
	par_buf.prec_gr[i] = itmp[2];
      }
    }

    /* Max restarts for cleanup solves */
    IF_OK par_buf.nrestart = 5;
    
    /* Precision for fermion force calculation */
    IF_OK status = get_i(stdin, prompt, "prec_ff", &par_buf.prec_ff);

    /*------------------------------------------------------------*/
    /* Chiral condensate and related quantities                   */
    /*------------------------------------------------------------*/
    
    IF_OK status += get_i(stdin, prompt, "number_of_pbp_masses",
			  &par_buf.num_pbp_masses);
    if(par_buf.num_pbp_masses > MAX_MASS_PBP){
      printf("Number of masses exceeds dimension %d\n",MAX_MASS_PBP);
      status++;
    }
    IF_OK if(par_buf.num_pbp_masses > 0){
      IF_OK status += get_i(stdin, prompt, "max_cg_prop",
			    &par_buf.qic_pbp[0].max);
      IF_OK status += get_i(stdin, prompt, "max_cg_prop_restarts",
			    &par_buf.qic_pbp[0].nrestart);
      IF_OK status += get_i(stdin, prompt, "npbp_reps", &par_buf.npbp_reps );
      IF_OK status += get_i(stdin, prompt, "prec_pbp", &par_buf.prec_pbp);
      IF_OK for(i = 0; i < par_buf.num_pbp_masses; i++){
	IF_OK status += get_f(stdin, prompt, "mass", &par_buf.ksp_pbp[i].mass);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	IF_OK status += get_f(stdin, prompt, "naik_term_epsilon", 
			      &par_buf.ksp_pbp[i].naik_term_epsilon);
#endif
	par_buf.qic_pbp[i].min = 0;
	par_buf.qic_pbp[i].start_flag = 0;
	par_buf.qic_pbp[i].nsrc = 1;
	par_buf.qic_pbp[i].max = par_buf.qic_pbp[0].max;
	par_buf.qic_pbp[i].nrestart = par_buf.qic_pbp[0].nrestart;
	par_buf.qic_pbp[i].prec = par_buf.prec_pbp;
	IF_OK status += get_f(stdin, prompt, "error_for_propagator", &par_buf.qic_pbp[i].resid);
	IF_OK status += get_f(stdin, prompt, "rel_error_for_propagator", &par_buf.qic_pbp[i].relresid );
      }
    }

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
					 par_buf.startfile );
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
				       par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
				 par_buf.stringLFN );
    
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  if( par_buf.stopflag != 0 )return par_buf.stopflag;
  
  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  propinterval = par_buf.propinterval;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  for(i = 0; i< n_pseudo; i++){
    niter_md[i] = par_buf.niter_md[i];
    niter_fa[i] = par_buf.niter_fa[i];
    niter_gr[i] = par_buf.niter_gr[i];

    rsqmin_md[i] = par_buf.rsqmin_md[i];
    rsqmin_fa[i] = par_buf.rsqmin_fa[i];
    rsqmin_gr[i] = par_buf.rsqmin_gr[i];

    prec_md[i] = par_buf.prec_md[i];
    prec_fa[i] = par_buf.prec_fa[i];
    prec_gr[i] = par_buf.prec_gr[i];
  }
  prec_ff = par_buf.prec_ff;
  rsqprop = par_buf.rsqprop;
  epsilon = par_buf.epsilon;
  n_pseudo = par_buf.n_pseudo;
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  int isite, idir;
  site *s;
  FORALLSITES(isite,s) {
    for( idir=XUP;idir<=TUP;idir++ ) {
      lattice[isite].on_step_Y[idir] = 0;
      lattice[isite].on_step_W[idir] = 0;
      lattice[isite].on_step_V[idir] = 0;
      /* zero out information from previous time step
         if fresh lattice; keep everything for continuation */
      if( startflag != CONTINUE ){
        lattice[isite].phase_Y_previous[idir] = 0.0;
        lattice[isite].phase_Y[idir] = 0.0;
        lattice[isite].Vdet[idir] = 0.0;
        clear_su3mat( &(lattice[isite].Wlink[idir]) );
        clear_su3mat( &(lattice[isite].Wlink_previous[idir]) );
      }
    }
  }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */

#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  /* Add PBP quantities to the eps_naik table of unique Naik epsilon
     coefficients .  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* Contribution from the chiral condensate epsilons */
  for(i = 0; i < par_buf.num_pbp_masses; i++){
    par_buf.ksp_pbp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, &n_naiks, 
		    par_buf.ksp_pbp[i].naik_term_epsilon);
  }

#endif

  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile );
  /* if a lattice was read in, put in KS phases and AP boundary condition */

  phases_in = OFF;
  rephase( ON );
  
  /* Copy gauge links from site structure to field */
  /* (This is a transitional step.  As we upgrade the code, we will 
     read/construct the lattice directly in the field and drop the
     site structure */

  /* Set options for fermion links */
#ifdef DM_DU0
  /* We want to calculate both the links and their u0 derivatives */
  fermion_links_want_du0(1);
#endif
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )& defined(DM_DEPS)
  /* We want to calculate both the links and their Naik eps
     derivatives (HISQ only) */
  fermion_links_want_deps(1);
#endif
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fn_links = create_fermion_links_from_site(PRECISION, n_naiks, eps_naik);
#else
  fn_links = create_fermion_links_from_site(PRECISION, 0, NULL);
#endif
    
  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();
  printf("Amit MyFFApp/setup.c end of function int readin(int prompt)\n");
  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
static void  make_3n_gathers(void)
{
  node0_printf("Amit MyFFApp/setup.c start of function static void make_3n_gathers \n");
  int i;
  
  for(i=XUP; i<=TUP; i++)
    {
      make_gather(third_neighbor, &i, WANT_INVERSE, ALLOW_EVEN_ODD, SWITCH_PARITY);
    }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */  
  sort_eight_gathers(X3UP);
  node0_printf("Amit MyFFApp/setup.c end of function static void make_3n_gathers \n");
  
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

static void  third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
	"forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  //node0_printf("Amit MyFFApp/setup.c start of function static void  third_neighbor(int x, int y, int z, int t, int *dirpt, int FB, int *xp, int *yp, int *zp, int *tp) \n");
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir)
    {
    case XUP: *xp = (x+3)%nx; break;
    case XDOWN: *xp = (x+4*nx-3)%nx; break;
    case YUP: *yp = (y+3)%ny; break;
    case YDOWN: *yp = (y+4*ny-3)%ny; break;
    case ZUP: *zp = (z+3)%nz; break;
    case ZDOWN: *zp = (z+4*nz-3)%nz; break;
    case TUP: *tp = (t+3)%nt; break;
    case TDOWN: *tp = (t+4*nt-3)%nt; break;
    default: printf("default MyFFApp/setup.c third_neighb: bad direction\n"); exit(1);
    }
  //node0_printf("Amit MyFFApp/setup.c end of function static void  third_neighbor(int x, int y, int z, int t, int *dirpt, int FB, int *xp, int *yp, int *zp, int *tp) \n");
}
