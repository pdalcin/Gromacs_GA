/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        single
 * Threads:          no
 * Software invsqrt: yes
 * PowerPC invsqrt:  no
 * Prefetch forces:  no
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<math.h>
#include<vec.h>



/*
 * Gromacs nonbonded kernel nb_kernel203
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Not calculated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel203(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         Vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float *         work)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         fscal,tx,ty,tz;
    float         rinvsq;
    float         jq;
    float         qq,vcoul,vctot;
    float         krsq;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         ix4,iy4,iz4,fix4,fiy4,fiz4;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx41,dy41,dz41,rsq41,rinv41;
    float         qH,qM;
    const int     fractshift = 12;
    const int     fractmask = 8388607;
    const int     expshift = 23;
    const int     expmask = 2139095040;
    const int     explsb = 8388608;
    float         lu;
    int           iexp,addr;
    union { unsigned int bval; float fval; } bitpattern,result;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];

    nj1              = 0;              
    
    for(n=0; (n<nri); n++)
    {
        is3              = 3*shift[n];     
        shX              = shiftvec[is3];  
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];      
        nj1              = jindex[n+1];    
        ii               = iinr[n];        
        ii3              = 3*ii;           
        ix2              = shX + pos[ii3+3];
        iy2              = shY + pos[ii3+4];
        iz2              = shZ + pos[ii3+5];
        ix3              = shX + pos[ii3+6];
        iy3              = shY + pos[ii3+7];
        iz3              = shZ + pos[ii3+8];
        ix4              = shX + pos[ii3+9];
        iy4              = shY + pos[ii3+10];
        iz4              = shZ + pos[ii3+11];
        vctot            = 0;              
        fix2             = 0;              
        fiy2             = 0;              
        fiz2             = 0;              
        fix3             = 0;              
        fiy3             = 0;              
        fiz3             = 0;              
        fix4             = 0;              
        fiy4             = 0;              
        fiz4             = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx21             = ix2 - jx1;      
            dy21             = iy2 - jy1;      
            dz21             = iz2 - jz1;      
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
            dx31             = ix3 - jx1;      
            dy31             = iy3 - jy1;      
            dz31             = iz3 - jz1;      
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
            dx41             = ix4 - jx1;      
            dy41             = iy4 - jy1;      
            dz41             = iz4 - jz1;      
            rsq41            = dx41*dx41+dy41*dy41+dz41*dz41;
            bitpattern.fval  = rsq21;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv21           = (0.5*lu*(3.0-((rsq21*lu)*lu)));
            bitpattern.fval  = rsq31;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv31           = (0.5*lu*(3.0-((rsq31*lu)*lu)));
            bitpattern.fval  = rsq41;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv41           = (0.5*lu*(3.0-((rsq41*lu)*lu)));
            jq               = charge[jnr+0];  
            qq               = qH*jq;          
            rinvsq           = rinv21*rinv21;  
            krsq             = krf*rsq21;      
            vcoul            = qq*(rinv21+krsq-crf);
            vctot            = vctot+vcoul;    
            fscal            = (qq*(rinv21-2.0*krsq))*rinvsq;
            tx               = fscal*dx21;     
            ty               = fscal*dy21;     
            tz               = fscal*dz21;     
            fix2             = fix2 + tx;      
            fiy2             = fiy2 + ty;      
            fiz2             = fiz2 + tz;      
            fjx1             = faction[j3+0] - tx;
            fjy1             = faction[j3+1] - ty;
            fjz1             = faction[j3+2] - tz;
            rinvsq           = rinv31*rinv31;  
            krsq             = krf*rsq31;      
            vcoul            = qq*(rinv31+krsq-crf);
            vctot            = vctot+vcoul;    
            fscal            = (qq*(rinv31-2.0*krsq))*rinvsq;
            tx               = fscal*dx31;     
            ty               = fscal*dy31;     
            tz               = fscal*dz31;     
            fix3             = fix3 + tx;      
            fiy3             = fiy3 + ty;      
            fiz3             = fiz3 + tz;      
            fjx1             = fjx1 - tx;      
            fjy1             = fjy1 - ty;      
            fjz1             = fjz1 - tz;      
            qq               = qM*jq;          
            rinvsq           = rinv41*rinv41;  
            krsq             = krf*rsq41;      
            vcoul            = qq*(rinv41+krsq-crf);
            vctot            = vctot+vcoul;    
            fscal            = (qq*(rinv41-2.0*krsq))*rinvsq;
            tx               = fscal*dx41;     
            ty               = fscal*dy41;     
            tz               = fscal*dz41;     
            fix4             = fix4 + tx;      
            fiy4             = fiy4 + ty;      
            fiz4             = fiz4 + tz;      
            faction[j3+0]    = fjx1 - tx;      
            faction[j3+1]    = fjy1 - ty;      
            faction[j3+2]    = fjz1 - tz;      
        }
        
        faction[ii3+3]   = faction[ii3+3] + fix2;
        faction[ii3+4]   = faction[ii3+4] + fiy2;
        faction[ii3+5]   = faction[ii3+5] + fiz2;
        faction[ii3+6]   = faction[ii3+6] + fix3;
        faction[ii3+7]   = faction[ii3+7] + fiy3;
        faction[ii3+8]   = faction[ii3+8] + fiz3;
        faction[ii3+9]   = faction[ii3+9] + fix4;
        faction[ii3+10]  = faction[ii3+10] + fiy4;
        faction[ii3+11]  = faction[ii3+11] + fiz4;
        fshift[is3]      = fshift[is3]+fix2+fix3+fix4;
        fshift[is3+1]    = fshift[is3+1]+fiy2+fiy3+fiy4;
        fshift[is3+2]    = fshift[is3+2]+fiz2+fiz3+fiz4;
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}





/*
 * Gromacs nonbonded kernel nb_kernel203nf
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Not calculated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        no
 */
void nb_kernel203nf(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         Vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float *         work)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         jq;
    float         qq,vcoul,vctot;
    float         krsq;
    float         ix2,iy2,iz2;
    float         ix3,iy3,iz3;
    float         ix4,iy4,iz4;
    float         jx1,jy1,jz1;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx41,dy41,dz41,rsq41,rinv41;
    float         qH,qM;
    const int     fractshift = 12;
    const int     fractmask = 8388607;
    const int     expshift = 23;
    const int     expmask = 2139095040;
    const int     explsb = 8388608;
    float         lu;
    int           iexp,addr;
    union { unsigned int bval; float fval; } bitpattern,result;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];

    nj1              = 0;              
    
    for(n=0; (n<nri); n++)
    {
        is3              = 3*shift[n];     
        shX              = shiftvec[is3];  
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];      
        nj1              = jindex[n+1];    
        ii               = iinr[n];        
        ii3              = 3*ii;           
        ix2              = shX + pos[ii3+3];
        iy2              = shY + pos[ii3+4];
        iz2              = shZ + pos[ii3+5];
        ix3              = shX + pos[ii3+6];
        iy3              = shY + pos[ii3+7];
        iz3              = shZ + pos[ii3+8];
        ix4              = shX + pos[ii3+9];
        iy4              = shY + pos[ii3+10];
        iz4              = shZ + pos[ii3+11];
        vctot            = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx21             = ix2 - jx1;      
            dy21             = iy2 - jy1;      
            dz21             = iz2 - jz1;      
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
            dx31             = ix3 - jx1;      
            dy31             = iy3 - jy1;      
            dz31             = iz3 - jz1;      
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
            dx41             = ix4 - jx1;      
            dy41             = iy4 - jy1;      
            dz41             = iz4 - jz1;      
            rsq41            = dx41*dx41+dy41*dy41+dz41*dz41;
            bitpattern.fval  = rsq21;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv21           = (0.5*lu*(3.0-((rsq21*lu)*lu)));
            bitpattern.fval  = rsq31;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv31           = (0.5*lu*(3.0-((rsq31*lu)*lu)));
            bitpattern.fval  = rsq41;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv41           = (0.5*lu*(3.0-((rsq41*lu)*lu)));
            jq               = charge[jnr+0];  
            qq               = qH*jq;          
            krsq             = krf*rsq21;      
            vcoul            = qq*(rinv21+krsq-crf);
            vctot            = vctot+vcoul;    
            krsq             = krf*rsq31;      
            vcoul            = qq*(rinv31+krsq-crf);
            vctot            = vctot+vcoul;    
            qq               = qM*jq;          
            krsq             = krf*rsq41;      
            vcoul            = qq*(rinv41+krsq-crf);
            vctot            = vctot+vcoul;    
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}


