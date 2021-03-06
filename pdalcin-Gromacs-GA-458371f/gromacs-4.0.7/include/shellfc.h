/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

#include "typedefs.h"
 
/* Initialization function, also predicts the initial shell postions.
 * If x!=NULL, the shells are predict for the global coordinates x.
 */
extern gmx_shellfc_t init_shell_flexcon(FILE *log,
					gmx_mtop_t *mtop,int nflexcon,
					rvec *x);

/* Get the local shell with domain decomposition */
extern void make_local_shells(t_commrec *cr,t_mdatoms *md,
			      gmx_shellfc_t shfc);

/* Optimize shell positions */
extern int relax_shell_flexcon(FILE *log,t_commrec *cr,bool bVerbose,
			       int mdstep,t_inputrec *inputrec,
			       bool bDoNS,bool bStopCM,
			       gmx_localtop_t *top,gmx_constr_t constr,
			       gmx_enerdata_t *enerd,t_fcdata *fcd,
			       t_state *state,rvec f[],
			       rvec buf[],tensor force_vir,
			       t_mdatoms *md,
			       t_nrnb *nrnb,gmx_wallcycle_t wcycle,
			       t_graph *graph,
			       gmx_groups_t *groups,
			       gmx_shellfc_t shfc,
			       t_forcerec *fr,
			       real t,rvec mu_tot,
			       int natoms,bool *bConverged,
			       gmx_vsite_t *vsite,
			       FILE *fp_field);
