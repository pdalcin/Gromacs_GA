/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "confio.h"
#include "copyrite.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "string2.h"
#include "vec.h"
#include "gstat.h"
#include "statutil.h"
#include "typedefs.h"
#include "xvgr.h"
	
static char *etitles[] = { "E-docked", "Free Energy" };

typedef struct {
  real    edocked,efree;
  int     index,cluster_id;
  t_atoms atoms;
  rvec    *x;
  int     ePBC;
  matrix  box;
} t_pdbfile;

static t_pdbfile *read_pdbf(char *fn)
{
  t_pdbfile *pdbf;
  double    e;
  char      buf[256],*ptr;
  int       natoms;
  FILE      *fp;
  
  snew(pdbf,1);
  get_stx_coordnum (fn,&natoms);
  init_t_atoms(&(pdbf->atoms),natoms,FALSE);
  snew(pdbf->x,natoms);
  read_stx_conf(fn,buf,&pdbf->atoms,pdbf->x,NULL,&pdbf->ePBC,pdbf->box);
  fp = fopen(fn,"r");
  do {
    ptr = fgets2(buf,255,fp);
    if (ptr) {
      if (strstr(buf,"Intermolecular") != NULL) {
	ptr = strchr(buf,'=');
	sscanf(ptr+1,"%lf",&e);
	pdbf->edocked = e;
      }
      else if (strstr(buf,"Estimated Free") != NULL) {
	ptr = strchr(buf,'=');
	sscanf(ptr+1,"%lf",&e);
	pdbf->efree = e;
      } 
    }
  } while (ptr != NULL);
  fclose(fp);
  
  return pdbf;
}

static t_pdbfile **read_em_all(char *fn,int *npdbf)
{
  t_pdbfile **pdbf=0;
  int  i,maxpdbf;
  char buf[256],name[256];
  bool bExist;
  
  strcpy(buf,fn);
  buf[strlen(buf)-4] = '\0';
  maxpdbf = 100;
  snew(pdbf,maxpdbf);
  i=0;
  do {
    sprintf(name,"%s_%d.pdb",buf,i+1);
    if ((bExist = fexist(name)) == TRUE) {
      pdbf[i] = read_pdbf(name);
      pdbf[i]->index = i+1;
      i++;
      if (i >= maxpdbf) {
	maxpdbf += 100;
	srenew(pdbf,maxpdbf);
      }
    }
  } while (bExist);
  
  *npdbf = i;

  printf("Found %d pdb files\n",i);
    
  return pdbf;
}

static bool bFreeSort=FALSE;

static int pdbf_comp(const void *a,const void *b)
{
  t_pdbfile *pa,*pb;
  real x;
  int  dc;
  
  pa = *(t_pdbfile **)a;
  pb = *(t_pdbfile **)b;
  
  dc = pa->cluster_id - pb->cluster_id;
  
  if (dc == 0) {
    if (bFreeSort)
      x = pa->efree   - pb->efree;
    else
      x = pa->edocked - pb->edocked;
    
    if (x < 0)
      return -1;
    else if (x > 0)
      return 1;
    else
      return 0;
  }
  else 
    return dc;
}

static void analyse_em_all(int npdb,t_pdbfile *pdbf[],
			   char *edocked,char *efree)
{
  FILE *fp;
  int i;
  
  for(bFreeSort = FALSE; (bFreeSort <= TRUE); bFreeSort++) {
    qsort(pdbf,npdb,sizeof(pdbf[0]),pdbf_comp);
    fp = xvgropen(bFreeSort ? efree : edocked,
		  etitles[bFreeSort],"()","E (kJ/mol)");
    for(i=0; (i<npdb); i++)
      fprintf(fp,"%12lf\n",bFreeSort ? pdbf[i]->efree : pdbf[i]->edocked);
    fclose(fp);
  }
}

static void clust_stat(FILE *fp,int start,int end,t_pdbfile *pdbf[])
{
  int i;
  t_lsq ed,ef;
  
  init_lsq(&ed);
  init_lsq(&ef);
  for(i=start; (i<end); i++) {
    add_lsq(&ed,i-start,pdbf[i]->edocked);
    add_lsq(&ef,i-start,pdbf[i]->efree);
  }
  fprintf(fp,"  <%12s> = %8.3f (+/- %6.3f)\n",etitles[FALSE],
	  aver_lsq(&ed),sigma_lsq(&ed));
  fprintf(fp,"  <%12s> = %8.3f (+/- %6.3f)\n",etitles[TRUE],
	  aver_lsq(&ef),sigma_lsq(&ef));
  done_lsq(&ed);
  done_lsq(&ef);
}

static real rmsd_dist(t_pdbfile *pa,t_pdbfile *pb,bool bRMSD)
{
  int  i;
  real rmsd,dist;
  rvec acm,bcm,dx;
  
  if (bRMSD) {
    rmsd = 0;
    for(i=0; (i<pa->atoms.nr); i++) {
      rvec_sub(pa->x[i],pb->x[i],dx);
      rmsd += iprod(dx,dx);
    }
    rmsd = sqrt(rmsd/pa->atoms.nr);
  }
  else {
    dist = 0;
    clear_rvec(acm);
    clear_rvec(bcm);
    for(i=0; (i<pa->atoms.nr); i++) {
      rvec_inc(acm,pa->x[i]);
      rvec_inc(bcm,pb->x[i]);
    }
    rvec_sub(acm,bcm,dx);
    for(i=0; (i<DIM); i++)
      dx[i] /= pa->atoms.nr;
    rmsd = norm(dx);
  }
  return rmsd;
}

static void line(FILE *fp)
{
  fprintf(fp,"                   - - - - - - - - - -\n");
}

static void cluster_em_all(FILE *fp,int npdb,t_pdbfile *pdbf[],
			   char *pdbout,bool bFree,bool bRMSD,real cutoff)
{
  int  i,j,k;
  int  *cndx,ncluster;
  real rmsd;
  
  bFreeSort = bFree;
  qsort(pdbf,npdb,sizeof(pdbf[0]),pdbf_comp);

  fprintf(fp,"Statistics over all structures:\n");
  clust_stat(fp,0,npdb,pdbf);
  line(fp);
  
  /* Index to first structure in a cluster */
  snew(cndx,npdb);
  ncluster=1;
  
  for(i=1; (i<npdb); i++) {
    for(j=0; (j<ncluster); j++) {
      rmsd = rmsd_dist(pdbf[cndx[j]],pdbf[i],bRMSD);
      if (rmsd <= cutoff) {
	/* Structure i is in cluster j */
	pdbf[i]->cluster_id = pdbf[cndx[j]]->cluster_id;
	break;
      }
    }
    if (j == ncluster) {
      /* New cluster! Cool! */
      cndx[ncluster] = i;
      pdbf[i]->cluster_id = ncluster;
      ncluster++;
    }
  }
  cndx[ncluster] = npdb;
  qsort(pdbf,npdb,sizeof(pdbf[0]),pdbf_comp);
  
  j=0;
  cndx[j++] = 0;
  for(i=1; (i<npdb); i++) {
    if (pdbf[i]->cluster_id != pdbf[i-1]->cluster_id) {
      cndx[j++] = i;
    }
  }
  cndx[j] = npdb;
  if (j != ncluster) 
    gmx_fatal(FARGS,"Consistency error: j = %d, ncluster = %d",j,ncluster);
  
  fprintf(fp,"I found %d clusters based on %s and %s with a %.3f nm cut-off\n",
	  ncluster,etitles[bFree],bRMSD ? "RMSD" : "distance",cutoff);
  line(fp);
  for(j=0; (j<ncluster); j++) {
    fprintf(fp,"Cluster: %3d  %s: %10.5f kJ/mol %3d elements\n",
	    j,etitles[bFree],
	    bFree ? pdbf[cndx[j]]->efree : pdbf[cndx[j]]->edocked,
	    cndx[j+1]-cndx[j]);
    clust_stat(fp,cndx[j],cndx[j+1],pdbf);
    for(k=cndx[j]; (k<cndx[j+1]); k++)
      fprintf(fp,"  %3d",pdbf[k]->index);
    fprintf(fp,"\n");
    line(fp);
  }
  sfree(cndx);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "anadock analyses the results of an Autodock run and clusters the",
    "structures together, based on distance or RMSD. The docked energy",
    "and free energy estimates are analysed, and for each cluster the",
    "energy statistics are printed.[PAR]",
    "An alternative approach to this is to cluster the structures first",
    "(using [TT]g_cluster[tt] and then sort the clusters on either lowest",
    "energy or average energy."
  };
  t_filenm fnm[] = {
    { efPDB, "-f", NULL,       ffREAD  },
    { efPDB, "-ox", "cluster", ffWRITE },
    { efXVG, "-od", "edocked", ffWRITE },
    { efXVG, "-of", "efree",   ffWRITE },
    { efLOG, "-g",  "anadock", ffWRITE }
  };
#define NFILE asize(fnm)
  static bool bFree=FALSE,bRMS=TRUE;
  static real cutoff = 0.2;
  t_pargs pa[] = {
    { "-free",   FALSE, etBOOL, {&bFree}, 
      "Use Free energy estimate from autodock for sorting the classes" },
    { "-rms",    FALSE, etBOOL, {&bRMS},
      "Cluster on RMS or distance" },
    { "-cutoff", FALSE, etREAL, {&cutoff},
      "Maximum RMSD/distance for belonging to the same cluster" }
  };
#define NPA asize(pa)

  FILE      *fp;
  t_pdbfile **pdbf=NULL;
  int       npdbf;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,NPA,pa,
		    asize(desc),desc,0,NULL);
  
  fp = ffopen(opt2fn("-g",NFILE,fnm),"w");
  please_cite(stdout,"Hetenyi2002b");
  please_cite(fp,"Hetenyi2002b");
  
  pdbf = read_em_all(opt2fn("-f",NFILE,fnm),&npdbf);
  
  analyse_em_all(npdbf,pdbf,opt2fn("-od",NFILE,fnm),opt2fn("-of",NFILE,fnm));
  
  cluster_em_all(fp,npdbf,pdbf,opt2fn("-ox",NFILE,fnm),
		 bFree,bRMS,cutoff);
  
  thanx(fp);
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
