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
 * Gromacs Runs On Most of All Computer Systems
 */

#include "mpelogging.h"

#ifndef _vec_h
#define _vec_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#define gmx_inline inline
#else
#ifdef __GNUC__
#define gmx_inline __inline
#else
#define inline
#endif
#endif
/*
  collection of in-line ready operations:
  
  lookup-table optimized scalar operations:
  real invsqrt(real x)
  void vecinvsqrt(real in[],real out[],int n)
  void vecrecip(real in[],real out[],int n)
  real sqr(real x)
  double dsqr(double x)
  
  vector operations:
  void rvec_add(const rvec a,const rvec b,rvec c)  c = a + b
  void dvec_add(const dvec a,const dvec b,dvec c)  c = a + b
  void ivec_add(const ivec a,const ivec b,ivec c)  c = a + b
  void rvec_inc(rvec a,const rvec b)               a += b
  void dvec_inc(dvec a,const dvec b)               a += b
  void ivec_inc(ivec a,const ivec b)               a += b
  void rvec_sub(const rvec a,const rvec b,rvec c)  c = a - b
  void dvec_sub(const dvec a,const dvec b,dvec c)  c = a - b
  void rvec_dec(rvec a,rvec b)                     a -= b
  void copy_rvec(const rvec a,rvec b)              b = a (reals)
  void copy_dvec(const dvec a,dvec b)              b = a (reals)
  void copy_ivec(const ivec a,ivec b)              b = a (integers)
  void ivec_sub(const ivec a,const ivec b,ivec c)  c = a - b
  void svmul(real a,rvec v1,rvec v2)               v2 = a * v1
  void dsvmul(double a,dvec v1,dvec v2)            v2 = a * v1
  void clear_rvec(rvec a)                          a = 0
  void clear_dvec(dvec a)                          a = 0
  void clear_ivec(rvec a)                          a = 0
  void clear_rvecs(int n,rvec v[])
  real iprod(rvec a,rvec b)                        = a . b (inner product)
  double diprod(dvec a,dvec b)                     = a . b (inner product)
  real iiprod(ivec a,ivec b)                       = a . b (integers)
  real norm2(rvec a)                               = | a |^2 ( = x*y*z )
  real norm(rvec a)                                = | a |
  double dnorm(dvec a)                             = | a |
  void cprod(rvec a,rvec b,rvec c)                 c = a x b (cross product)
  void dprod(rvec a,rvec b,rvec c)                 c = a * b (direct product)
  real cos_angle(rvec a,rvec b)
  real cos_angle_no_table(rvec a,rvec b)
  real distance2(rvec v1, rvec v2)                 = | v2 - v1 |^2
  void unitv(rvec src,rvec dest)                   dest = src / |src|
  void unitv_no_table(rvec src,rvec dest)          dest = src / |src|
  
  matrix (3x3) operations:
    ! indicates that dest should not be the same as a, b or src
    the _ur0 varieties work on matrices that have only zeros
    in the upper right part, such as box matrices, these varieties
    could produce less rounding errors, not due to the operations themselves,
    but because the compiler can easier recombine the operations
  void copy_mat(matrix a,matrix b)                 b = a
  void clear_mat(matrix a)			   a = 0
  void mmul(matrix a,matrix b,matrix dest)	!  dest = a . b
  void mmul_ur0(matrix a,matrix b,matrix dest)     dest = a . b
  void transpose(matrix src,matrix dest)	!  dest = src*
  void tmmul(matrix a,matrix b,matrix dest)	!  dest = a* . b
  void mtmul(matrix a,matrix b,matrix dest)	!  dest = a . b*
  real det(matrix a)				   = det(a)
  void m_add(matrix a,matrix b,matrix dest)	   dest = a + b
  void m_sub(matrix a,matrix b,matrix dest)	   dest = a - b
  void msmul(matrix m1,real r1,matrix dest)	   dest = r1 * m1
  void m_inv_ur0(matrix src,matrix dest)           dest = src^-1
  void m_inv(matrix src,matrix dest)		!  dest = src^-1
  void mvmul(matrix a,rvec src,rvec dest)	!  dest = a . src
  void mvmul_ur0(matrix a,rvec src,rvec dest)      dest = a . src
  void tmvmul_ur0(matrix a,rvec src,rvec dest)     dest = a* . src
  real trace(matrix m)                             = trace(m)
*/

#include "maths.h"
#include "typedefs.h"
#include "sysstuff.h"
#include "macros.h"
#include "gmx_fatal.h"

#define EXP_LSB         0x00800000
#define EXP_MASK        0x7f800000
#define EXP_SHIFT       23
#define FRACT_MASK      0x007fffff
#define FRACT_SIZE      11              /* significant part of fraction */
#define FRACT_SHIFT     (EXP_SHIFT-FRACT_SIZE)
#define EXP_ADDR(val)   (((val)&EXP_MASK)>>EXP_SHIFT)
#define FRACT_ADDR(val) (((val)&(FRACT_MASK|EXP_LSB))>>FRACT_SHIFT)

#define PR_VEC(a)       a[XX],a[YY],a[ZZ]

#ifdef GMX_SOFTWARE_SQRT
extern const unsigned int *  gmx_invsqrt_exptab;
extern const unsigned int *  gmx_invsqrt_fracttab;
#endif


typedef union 
{
  unsigned int bval;
  float fval;
} t_convert;


#ifdef GMX_SOFTWARE_SQRT
static inline real invsqrt(real x)
{
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  real        lu;
  real        y;
#ifdef GMX_DOUBLE
  real        y2;
#endif
 
  bit_pattern.fval=x;
  exp   = EXP_ADDR(bit_pattern.bval);
  fract = FRACT_ADDR(bit_pattern.bval);
  result.bval=gmx_invsqrt_exptab[exp] | gmx_invsqrt_fracttab[fract];
  lu    = result.fval;
  
  y=(half*lu*(three-((x*lu)*lu)));
#ifdef GMX_DOUBLE
  y2=(half*y*(three-((x*y)*y)));
  
  return y2;                    /* 10 Flops */
#else
  return y;                     /* 5  Flops */
#endif
}
#define INVSQRT_DONE 
#endif /* gmx_invsqrt */

#ifdef GMX_POWERPC_SQRT
static inline real invsqrt(real x)
{
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  real        lu;
  real        y;
#ifdef GMX_DOUBLE
  real        y2;
#endif

  lu = __frsqrte((double)x);

  y=(half*lu*(three-((x*lu)*lu)));

#if (GMX_POWERPC_SQRT==2)
  /* Extra iteration required */
  y=(half*y*(three-((x*y)*y)));
#endif

#ifdef GMX_DOUBLE
  y2=(half*y*(three-((x*y)*y)));

  return y2;                    /* 10 Flops */
#else
  return y;                     /* 5  Flops */
#endif
}
#define INVSQRT_DONE
#endif /* powerpc_invsqrt */


#ifndef INVSQRT_DONE
#define invsqrt(x) (1.0f/sqrt(x))
#endif





static inline real sqr(real x)
{
  return (x*x);
}

static inline double dsqr(double x)
{
  return (x*x);
}

extern void vecinvsqrt(real in[],real out[],int n);
/* Perform out[i]=1.0/sqrt(in[i]) for n elements */


extern void vecrecip(real in[],real out[],int n);
/* Perform out[i]=1.0/(in[i]) for n elements */

/* Note: If you need a fast version of vecinvsqrt 
 * and/or vecrecip, call detectcpu() and run the SSE/3DNow/SSE2/Altivec
 * versions if your hardware supports it.
 *
 * To use those routines, your memory HAS TO BE CACHE-ALIGNED.
 * Start by allocating 31 bytes more than you need, put
 * this in a temp variable (e.g. _buf, so you can free it later), and
 * create your aligned array buf with
 * 
 *  buf=(real *) ( ( (unsigned long int)_buf + 31 ) & (~0x1f) );
 */


static inline void rvec_add(const rvec a,const rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void dvec_add(const dvec a,const dvec b,dvec c)
{
  double x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void ivec_add(const ivec a,const ivec b,ivec c)
{
  int x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void rvec_inc(rvec a,const rvec b)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static inline void dvec_inc(dvec a,const dvec b)
{
  double x,y,z;

  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];

  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static inline void rvec_sub(const rvec a,const rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void dvec_sub(const dvec a,const dvec b,dvec c)
{
  double x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void rvec_dec(rvec a,const rvec b)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static inline void copy_rvec(const rvec a,rvec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static inline void copy_dvec(const dvec a,dvec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static inline void copy_ivec(const ivec a,ivec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static inline void ivec_sub(const ivec a,const ivec b,ivec c)
{
  int x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void copy_mat(matrix a,matrix b)
{
  copy_rvec(a[XX],b[XX]);
  copy_rvec(a[YY],b[YY]);
  copy_rvec(a[ZZ],b[ZZ]);
}

static inline void svmul(real a,const rvec v1,rvec v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
}

static inline void dsvmul(double a,const dvec v1,dvec v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
}

static inline real distance2(const rvec v1,const rvec v2)
{
  return sqr(v2[XX]-v1[XX]) + sqr(v2[YY]-v1[YY]) + sqr(v2[ZZ]-v1[ZZ]);
}

static inline void clear_rvec(rvec a)
{
  /* The ibm compiler has problems with inlining this 
   * when we use a const real variable
   */
  a[XX]=0.0;
  a[YY]=0.0;
  a[ZZ]=0.0;
}

static inline void clear_dvec(dvec a)
{
  /* The ibm compiler has problems with inlining this 
   * when we use a const real variable
   */
  a[XX]=0.0;
  a[YY]=0.0;
  a[ZZ]=0.0;
}

static inline void clear_ivec(ivec a)
{
  a[XX]=0;
  a[YY]=0;
  a[ZZ]=0;
}

static inline void clear_rvecs(int n,rvec v[])
{
/*  memset(v[0],0,DIM*n*sizeof(v[0][0])); */
  int i;

  GMX_MPE_LOG(ev_clear_rvecs_start);
    
  for(i=0; (i<n); i++) 
    clear_rvec(v[i]);
    
  GMX_MPE_LOG(ev_clear_rvecs_finish);  
}

static inline void clear_mat(matrix a)
{
/*  memset(a[0],0,DIM*DIM*sizeof(a[0][0])); */
  
  const real nul=0.0;
  
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=nul;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=nul;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=nul;
}

static inline real iprod(const rvec a,const rvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline double diprod(const dvec a,const dvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real iiprod(const ivec a,const ivec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real norm2(const rvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static inline real norm(const rvec a)
{
  return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

static inline double dnorm(const dvec a)
{
  return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

static inline real cos_angle(const rvec a,const rvec b)
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  real   cos;
  int    m;
  double aa,bb,ip,ipa,ipb,ipab; /* For accuracy these must be double! */
  
  ip=ipa=ipb=0.0;
  for(m=0; (m<DIM); m++) {		/* 18		*/
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  ipab = ipa*ipb;
  if (ipab > 0)
    cos = ip*invsqrt(ipab);		/*  7		*/
  else 
    cos = 1;
					/* 25 TOTAL	*/
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}

static inline real cos_angle_no_table(const rvec a,const rvec b)
{
  /* This version does not need the invsqrt lookup table */
  real   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */
  
  ip=ipa=ipb=0.0;
  for(m=0; (m<DIM); m++) {		/* 18		*/
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb); 		/* 12		*/
					/* 30 TOTAL	*/
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}

static inline void cprod(const rvec a,const rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static inline void mmul_ur0(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]*b[XX][XX];
  dest[XX][YY]=0.0;
  dest[XX][ZZ]=0.0;
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX];
  dest[YY][YY]=                    a[YY][YY]*b[YY][YY];
  dest[YY][ZZ]=0.0;
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[ZZ][YY]=                    a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[ZZ][ZZ]=                                        a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void mmul(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[XX][YY]*b[YY][XX]+a[XX][ZZ]*b[ZZ][XX];
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[YY][ZZ]*b[ZZ][XX];
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[XX][YY]=a[XX][XX]*b[XX][YY]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[ZZ][YY];
  dest[YY][YY]=a[YY][XX]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[ZZ][YY];
  dest[ZZ][YY]=a[ZZ][XX]*b[XX][YY]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[XX][ZZ]=a[XX][XX]*b[XX][ZZ]+a[XX][YY]*b[YY][ZZ]+a[XX][ZZ]*b[ZZ][ZZ];
  dest[YY][ZZ]=a[YY][XX]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[YY][ZZ]*b[ZZ][ZZ];
  dest[ZZ][ZZ]=a[ZZ][XX]*b[XX][ZZ]+a[ZZ][YY]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void transpose(matrix src,matrix dest)
{
  dest[XX][XX]=src[XX][XX];
  dest[YY][XX]=src[XX][YY];
  dest[ZZ][XX]=src[XX][ZZ];
  dest[XX][YY]=src[YY][XX];
  dest[YY][YY]=src[YY][YY];
  dest[ZZ][YY]=src[YY][ZZ];
  dest[XX][ZZ]=src[ZZ][XX];
  dest[YY][ZZ]=src[ZZ][YY];
  dest[ZZ][ZZ]=src[ZZ][ZZ];
}

static inline void tmmul(matrix a,matrix b,matrix dest)
{
  /* Computes dest=mmul(transpose(a),b,dest) - used in do_pr_pcoupl */
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[YY][XX]*b[YY][XX]+a[ZZ][XX]*b[ZZ][XX];
  dest[XX][YY]=a[XX][XX]*b[XX][YY]+a[YY][XX]*b[YY][YY]+a[ZZ][XX]*b[ZZ][YY];
  dest[XX][ZZ]=a[XX][XX]*b[XX][ZZ]+a[YY][XX]*b[YY][ZZ]+a[ZZ][XX]*b[ZZ][ZZ];
  dest[YY][XX]=a[XX][YY]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[ZZ][YY]*b[ZZ][XX];
  dest[YY][YY]=a[XX][YY]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[ZZ][YY]*b[ZZ][YY];
  dest[YY][ZZ]=a[XX][YY]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[ZZ][YY]*b[ZZ][ZZ];
  dest[ZZ][XX]=a[XX][ZZ]*b[XX][XX]+a[YY][ZZ]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[ZZ][YY]=a[XX][ZZ]*b[XX][YY]+a[YY][ZZ]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[ZZ][ZZ]=a[XX][ZZ]*b[XX][ZZ]+a[YY][ZZ]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void mtmul(matrix a,matrix b,matrix dest)
{
  /* Computes dest=mmul(a,transpose(b),dest) - used in do_pr_pcoupl */
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[XX][YY]*b[XX][YY]+a[XX][ZZ]*b[XX][ZZ];
  dest[XX][YY]=a[XX][XX]*b[YY][XX]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[YY][ZZ];
  dest[XX][ZZ]=a[XX][XX]*b[ZZ][XX]+a[XX][YY]*b[ZZ][YY]+a[XX][ZZ]*b[ZZ][ZZ];
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[XX][YY]+a[YY][ZZ]*b[XX][ZZ];
  dest[YY][YY]=a[YY][XX]*b[YY][XX]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[YY][ZZ];
  dest[YY][ZZ]=a[YY][XX]*b[ZZ][XX]+a[YY][YY]*b[ZZ][YY]+a[YY][ZZ]*b[ZZ][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[XX][YY]+a[ZZ][ZZ]*b[XX][ZZ];
  dest[ZZ][YY]=a[ZZ][XX]*b[YY][XX]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[YY][ZZ];
  dest[ZZ][ZZ]=a[ZZ][XX]*b[ZZ][XX]+a[ZZ][YY]*b[ZZ][YY]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline real det(matrix a)
{
  return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
	  -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
	  +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}

static inline void m_add(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]+b[XX][XX];
  dest[XX][YY]=a[XX][YY]+b[XX][YY];
  dest[XX][ZZ]=a[XX][ZZ]+b[XX][ZZ];
  dest[YY][XX]=a[YY][XX]+b[YY][XX];
  dest[YY][YY]=a[YY][YY]+b[YY][YY];
  dest[YY][ZZ]=a[YY][ZZ]+b[YY][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]+b[ZZ][XX];
  dest[ZZ][YY]=a[ZZ][YY]+b[ZZ][YY];
  dest[ZZ][ZZ]=a[ZZ][ZZ]+b[ZZ][ZZ];
}

static inline void m_sub(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]-b[XX][XX];
  dest[XX][YY]=a[XX][YY]-b[XX][YY];
  dest[XX][ZZ]=a[XX][ZZ]-b[XX][ZZ];
  dest[YY][XX]=a[YY][XX]-b[YY][XX];
  dest[YY][YY]=a[YY][YY]-b[YY][YY];
  dest[YY][ZZ]=a[YY][ZZ]-b[YY][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]-b[ZZ][XX];
  dest[ZZ][YY]=a[ZZ][YY]-b[ZZ][YY];
  dest[ZZ][ZZ]=a[ZZ][ZZ]-b[ZZ][ZZ];
}

static inline void msmul(matrix m1,real r1,matrix dest)
{
  dest[XX][XX]=r1*m1[XX][XX];
  dest[XX][YY]=r1*m1[XX][YY];
  dest[XX][ZZ]=r1*m1[XX][ZZ];
  dest[YY][XX]=r1*m1[YY][XX];
  dest[YY][YY]=r1*m1[YY][YY];
  dest[YY][ZZ]=r1*m1[YY][ZZ];
  dest[ZZ][XX]=r1*m1[ZZ][XX];
  dest[ZZ][YY]=r1*m1[ZZ][YY];
  dest[ZZ][ZZ]=r1*m1[ZZ][ZZ];
}

static inline void m_inv_ur0(matrix src,matrix dest)
{
  double tmp = src[XX][XX]*src[YY][YY]*src[ZZ][ZZ];
  if (gmx_within_tol(tmp,0.0,100*GMX_REAL_MIN))
    gmx_fatal(FARGS,"Can not invert matrix, determinant is zero");

  dest[XX][XX] = 1/src[XX][XX];
  dest[YY][YY] = 1/src[YY][YY];
  dest[ZZ][ZZ] = 1/src[ZZ][ZZ];
  dest[ZZ][XX] = (src[YY][XX]*src[ZZ][YY]*dest[YY][YY]
		  - src[ZZ][XX])*dest[XX][XX]*dest[ZZ][ZZ];
  dest[YY][XX] = -src[YY][XX]*dest[XX][XX]*dest[YY][YY];
  dest[ZZ][YY] = -src[ZZ][YY]*dest[YY][YY]*dest[ZZ][ZZ];
  dest[XX][YY] = 0.0;
  dest[XX][ZZ] = 0.0;
  dest[YY][ZZ] = 0.0;
}

static inline void m_inv(matrix src,matrix dest)
{
  const real smallreal = 1.0e-24;
  const real largereal = 1.0e24;
  real  deter,c,fc;

  deter = det(src);
  c     = 1.0/deter;
  fc    = fabs(c);
  
  if ((fc <= smallreal) || (fc >= largereal)) 
    gmx_fatal(FARGS,"Can not invert matrix, determinant = %e",deter);

  dest[XX][XX]= c*(src[YY][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[YY][ZZ]);
  dest[XX][YY]=-c*(src[XX][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[XX][ZZ]);
  dest[XX][ZZ]= c*(src[XX][YY]*src[YY][ZZ]-src[YY][YY]*src[XX][ZZ]);
  dest[YY][XX]=-c*(src[YY][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[YY][ZZ]);
  dest[YY][YY]= c*(src[XX][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[XX][ZZ]);
  dest[YY][ZZ]=-c*(src[XX][XX]*src[YY][ZZ]-src[YY][XX]*src[XX][ZZ]);
  dest[ZZ][XX]= c*(src[YY][XX]*src[ZZ][YY]-src[ZZ][XX]*src[YY][YY]);
  dest[ZZ][YY]=-c*(src[XX][XX]*src[ZZ][YY]-src[ZZ][XX]*src[XX][YY]);
  dest[ZZ][ZZ]= c*(src[XX][XX]*src[YY][YY]-src[YY][XX]*src[XX][YY]);
}

static inline void mvmul(matrix a,const rvec src,rvec dest)
{
  dest[XX]=a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
  dest[YY]=a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
  dest[ZZ]=a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}

static inline void mvmul_ur0(matrix a,const rvec src,rvec dest)
{
  dest[ZZ]=a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
  dest[YY]=a[YY][XX]*src[XX]+a[YY][YY];
  dest[XX]=a[XX][XX]*src[XX];
}

static inline void tmvmul_ur0(matrix a,const rvec src,rvec dest)
{
  dest[XX]=a[XX][XX]*src[XX]+a[YY][XX]*src[YY]+a[ZZ][XX]*src[ZZ];
  dest[YY]=                  a[YY][YY]*src[YY]+a[ZZ][YY]*src[ZZ];
  dest[ZZ]=                                    a[ZZ][ZZ]*src[ZZ];
}

static inline void unitv(const rvec src,rvec dest)
{
  real linv;
  
  linv=invsqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static inline void unitv_no_table(const rvec src,rvec dest)
{
  real linv;
  
  linv=1.0/sqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static void calc_lll(rvec box,rvec lll)
{
  lll[XX] = 2.0*M_PI/box[XX];
  lll[YY] = 2.0*M_PI/box[YY];
  lll[ZZ] = 2.0*M_PI/box[ZZ];
}

static inline real trace(matrix m)
{
  return (m[XX][XX]+m[YY][YY]+m[ZZ][ZZ]);
}

static inline real _divide(real a,real b,char *file,int line)
{
    if (gmx_within_tol(b,0.0,GMX_REAL_MIN)) 
        gmx_fatal(FARGS,"Dividing by zero, file %s, line %d",file,line);
    return a/b;
}

static inline int _mod(int a,int b,char *file,int line)
{
  if(b==0)
    gmx_fatal(FARGS,"Modulo zero, file %s, line %d",file,line);
  return a % b;
}

/* Math on multidimensional rvecs, used e.g. in edsam.c */
static void m_rvecsub(int dim, rvec *a, rvec *b, rvec *c) 
{
    /* c = a-b */
    int i;

    for (i=0; i<dim; i++)
        rvec_sub(a[i], b[i], c[i]);
}


static void m_rvecadd(int dim, rvec *a, rvec *b, rvec *c)
{
    /* c = a+b */
    int i;

    for (i=0; i<dim; i++)
        rvec_add(a[i], b[i], c[i]);
}


static real m_rvecnorm(int dim, rvec *a)
{
    /* c = a+b */
    int i;
    real sum=0;

    for (i=0; i<dim; i++)
        sum += norm2(a[i]);

    return sqrt(sum);
}


static void m_rvecsmul(int dim, real s, rvec *a)
{
    int i;

    for (i=0; i<dim; i++)
        svmul(s, a[i], a[i]);
} 


static void m_rvec_to_one(int dim, rvec *a)
{
    m_rvecsmul(dim, 1.0/m_rvecnorm(dim,a), a);
}


static void m_rveccopy(int dim, rvec *a, rvec *b)
{
    /* b = a */
    int i;

    for (i=0; i<dim; i++)
        copy_rvec(a[i],b[i]);
} 


#define divide(a,b) _divide((a),(b),__FILE__,__LINE__)
#define mod(a,b)    _mod((a),(b),__FILE__,__LINE__)

#endif	/* _vec_h */
