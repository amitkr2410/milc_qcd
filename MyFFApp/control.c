/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"

#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int main( int argc, char **argv )
{
  int i,MeasurementCount,traj_done, naik_index;
  int prompt;
  int s_iters, avs_iters, avbcorr_iters;
  double dtime, dclock();

  //Plaquette and Field-strength variable
  double ss_plaq=0.0, st_plaq=0.0;
  complex TraceF3iF3iMinusF4iF4i, TraceF4iF3iPlusF3iF4i;
  
  // Initialization 
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1); 
  g_sync();
  
  /* set up lattice parameters */
  prompt = setup();

  printf("Amit MyFFApp/control.c prompt= %d \n",prompt);
  printf("Amit MyFFApp/control.c before while(readin(prompt)==0) called \n");
  /* loop over input sets */
  while( readin(prompt) == 0)
    {
      
      /* perform warmup trajectories */
      #ifdef MILC_GLOBAL_DEBUG
      global_current_time_step = 0;
      #endif /* MILC_GLOBAL_DEBUG */
    
      // Start application timer
      dtime = -dclock();
      printf(" Amit MyFFApp/control.c inside while(readin(prompt)==0) \n");
      for( traj_done=0; traj_done < warms; traj_done++ )
      	{
	  //rephase(OFF);
	  ss_plaq=0.0; st_plaq=0.0;
	  d_plaquette(&ss_plaq, &st_plaq);
	  printf("Amit MyFFApp/control.c Plaquette = (%e,%e)\n",ss_plaq, st_plaq);	  
	  //rephase(ON);
      	  update();	  
        }
      
      node0_printf("default MyFFApp/control.c WARMUPS COMPLETED\n"); fflush(stdout);
      
      /* perform measuring trajectories, reunitarizing and measuring 	*/
      MeasurementCount=0;		/* number of measurements 		*/
      avs_iters = avbcorr_iters = 0;

      for( traj_done=0; traj_done < trajecs; traj_done++ )
	{ 
          #ifdef MILC_GLOBAL_DEBUG
          #ifdef HISQ_REUNITARIZATION_DEBUG
	  {
	    int isite, idir;
	    site *s;
	    FORALLSITES(isite,s)
	      {
		for( idir=XUP;idir<=TUP;idir++ )
		  {
		    lattice[isite].on_step_Y[idir] = 0;
		    lattice[isite].on_step_W[idir] = 0;
		    lattice[isite].on_step_V[idir] = 0;
		  }
	      }
	  }
          #endif /* HISQ_REUNITARIZATION_DEBUG */
          #endif /* MILC_GLOBAL_DEBUG */
	  printf(" Amit MyFFApp/control.c s_iters=update() will call \n");
	  
	  /* do the trajectories */
	  //rephase(ON);
	  s_iters=update();
          printf(" Amit MyFFApp/control.c s_iters=update() called \n");
	  /* measure every "propinterval" trajectories */
	 if( (traj_done%propinterval)==(propinterval-1) )
	  {	      
	      /* call gauge_variable fermion_variable measuring routines */
	      //rephase(OFF);	      
	      /* Compute plaquette and display output */
	      ss_plaq=0.0; st_plaq=0.0;
	      d_plaquette(&ss_plaq, &st_plaq);
	      printf("Amit MyFFApp/control.c Plaquette = (%e,%e)\n",ss_plaq, st_plaq);
	      
	      /* Calculate trace of fmunu and output */
	      TraceF3iF3iMinusF4iF4i.real = 0.0; TraceF3iF3iMinusF4iF4i.imag = 0.0;
	      TraceF4iF3iPlusF3iF4i.real = 0.0;  TraceF4iF3iPlusF3iF4i.imag = 0.0;
	      //rephase(OFF);
	      fmunu_fmunu(TraceF3iF3iMinusF4iF4i, TraceF4iF3iPlusF3iF4i);
	      printf("Amit MyFFApp/control.c TraceF3iF3iMinusF4iF4i = (%e,%e)\n",TraceF3iF3iMinusF4iF4i.real, TraceF3iF3iMinusF4iF4i.imag);
	      printf("Amit MyFFApp/control.c TraceF4iF3iPlusF3iF4i  = (%e,%e)\n",TraceF4iF3iPlusF3iF4i.real, TraceF4iF3iPlusF3iF4i.imag);
	      
	      avs_iters += s_iters;
	      MeasurementCount = MeasurementCount + 1;
	      fflush(stdout);
	     }
	}	/* end loop over trajectories */
    
      node0_printf("default MyFFApp/control.c RUNNING COMPLETED\n"); fflush(stdout);
      if(MeasurementCount>0)  {
	node0_printf("default MyFFApp/control.c average cg iters for step= %e\n",
		     (double)avs_iters/MeasurementCount);
      }
      
      dtime += dclock();
      if(this_node==0){
	printf("Default MyFFApp/control.c Time = %e seconds\n",dtime);
	printf("Default MyFFApp/control.c total_iters = %d\n",total_iters);
      }
      fflush(stdout);
      
      /* save lattice if requested */
      if( saveflag != FORGET ){
	rephase( OFF );
	save_lattice( saveflag, savefile, stringLFN );
	rephase( ON );
      }
      
      /* Destroy fermion links (created in readin() */
      
#if FERM_ACTION == HISQ
      destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
      destroy_fermion_links_hypisq(fn_links);
#else
      destroy_fermion_links(fn_links);
#endif
      fn_links = NULL;
    }
  
  
  normal_exit(0);
  return 0;
}

