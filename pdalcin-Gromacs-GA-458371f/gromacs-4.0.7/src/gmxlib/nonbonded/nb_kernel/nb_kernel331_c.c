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
 * Gromacs nonbonded kernel nb_kernel331
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Tabulated
 * water optimization:      SPC/TIP3P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel331(
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
    float         jq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         FF;
    float         fijC;
    float         fijD,fijR;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         qO,qH;
    float         c6,c12;
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
    qO               = facel*charge[ii];
    qH               = facel*charge[ii+1];
    nti              = 2*ntype*type[ii];

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
        ix1              = shX + pos[ii3+0];
        iy1              = shY + pos[ii3+1];
        iz1              = shZ + pos[ii3+2];
        ix2              = shX + pos[ii3+3];
        iy2              = shY + pos[ii3+4];
        iz2              = shZ + pos[ii3+5];
        ix3              = shX + pos[ii3+6];
        iy3              = shY + pos[ii3+7];
        iz3              = shZ + pos[ii3+8];
        vctot            = 0;              
        Vvdwtot          = 0;              
        fix1             = 0;              
        fiy1             = 0;              
        fiz1             = 0;              
        fix2             = 0;              
        fiy2             = 0;              
        fiz2             = 0;              
        fix3             = 0;              
        fiy3             = 0;              
        fiz3             = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
            dx21             = ix2 - jx1;      
            dy21             = iy2 - jy1;      
            dz21             = iz2 - jz1;      
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
            dx31             = ix3 - jx1;      
            dy31             = iy3 - jy1;      
            dz31             = iz3 - jz1;      
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
            bitpattern.fval  = rsq11;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv11           = (0.5*lu*(3.0-((rsq11*lu)*lu)));
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
            jq               = charge[jnr+0];  
            qq               = qO*jq;          
            tj               = nti+2*type[jnr];
            c6               = vdwparam[tj];   
            c12              = vdwparam[tj+1]; 
            r                = rsq11*rinv11;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            vcoul            = qq*VV;          
            fijC             = qq*FF;          
            vctot            = vctot + vcoul;  
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            Vvdw6            = c6*VV;          
            fijD             = c6*FF;          
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            Vvdw12           = c12*VV;         
            fijR             = c12*FF;         
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;
            fscal            = -((fijC+fijD+fijR)*tabscale)*rinv11;
            tx               = fscal*dx11;     
            ty               = fscal*dy11;     
            tz               = fscal*dz11;     
            fix1             = fix1 + tx;      
            fiy1             = fiy1 + ty;      
            fiz1             = fiz1 + tz;      
            fjx1             = faction[j3+0] - tx;
            fjy1             = faction[j3+1] - ty;
            fjz1             = faction[j3+2] - tz;
            qq               = qH*jq;          
            r                = rsq21*rinv21;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            vcoul            = qq*VV;          
            fijC             = qq*FF;          
            vctot            = vctot + vcoul;  
            fscal            = -((fijC)*tabscale)*rinv21;
            tx               = fscal*dx21;     
            ty               = fscal*dy21;     
            tz               = fscal*dz21;     
            fix2             = fix2 + tx;      
            fiy2             = fiy2 + ty;      
            fiz2             = fiz2 + tz;      
            fjx1             = fjx1 - tx;      
            fjy1             = fjy1 - ty;      
            fjz1             = fjz1 - tz;      
            r                = rsq31*rinv31;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            vcoul            = qq*VV;          
            fijC             = qq*FF;          
            vctot            = vctot + vcoul;  
            fscal            = -((fijC)*tabscale)*rinv31;
            tx               = fscal*dx31;     
            ty               = fscal*dy31;     
            tz               = fscal*dz31;     
            fix3             = fix3 + tx;      
            fiy3             = fiy3 + ty;      
            fiz3             = fiz3 + tz;      
            faction[j3+0]    = fjx1 - tx;      
            faction[j3+1]    = fjy1 - ty;      
            faction[j3+2]    = fjz1 - tz;      
        }
        
        faction[ii3+0]   = faction[ii3+0] + fix1;
        faction[ii3+1]   = faction[ii3+1] + fiy1;
        faction[ii3+2]   = faction[ii3+2] + fiz1;
        faction[ii3+3]   = faction[ii3+3] + fix2;
        faction[ii3+4]   = faction[ii3+4] + fiy2;
        faction[ii3+5]   = faction[ii3+5] + fiz2;
        faction[ii3+6]   = faction[ii3+6] + fix3;
        faction[ii3+7]   = faction[ii3+7] + fiy3;
        faction[ii3+8]   = faction[ii3+8] + fiz3;
        fshift[is3]      = fshift[is3]+fix1+fix2+fix3;
        fshift[is3+1]    = fshift[is3+1]+fiy1+fiy2+fiy3;
        fshift[is3+2]    = fshift[is3+2]+fiz1+fiz2+fiz3;
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}





/*
 * Gromacs nonbonded kernel nb_kernel331nf
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Tabulated
 * water optimization:      SPC/TIP3P - other atoms
 * Calculate forces:        no
 */
void nb_kernel331nf(
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
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         ix1,iy1,iz1;
    float         ix2,iy2,iz2;
    float         ix3,iy3,iz3;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         qO,qH;
    float         c6,c12;
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
    qO               = facel*charge[ii];
    qH               = facel*charge[ii+1];
    nti              = 2*ntype*type[ii];

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
        ix1              = shX + pos[ii3+0];
        iy1              = shY + pos[ii3+1];
        iz1              = shZ + pos[ii3+2];
        ix2              = shX + pos[ii3+3];
        iy2              = shY + pos[ii3+4];
        iz2              = shZ + pos[ii3+5];
        ix3              = shX + pos[ii3+6];
        iy3              = shY + pos[ii3+7];
        iz3              = shZ + pos[ii3+8];
        vctot            = 0;              
        Vvdwtot          = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
            dx21             = ix2 - jx1;      
            dy21             = iy2 - jy1;      
            dz21             = iz2 - jz1;      
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
            dx31             = ix3 - jx1;      
            dy31             = iy3 - jy1;      
            dz31             = iz3 - jz1;      
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
            bitpattern.fval  = rsq11;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv11           = (0.5*lu*(3.0-((rsq11*lu)*lu)));
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
            jq               = charge[jnr+0];  
            qq               = qO*jq;          
            tj               = nti+2*type[jnr];
            c6               = vdwparam[tj];   
            c12              = vdwparam[tj+1]; 
            r                = rsq11*rinv11;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            Vvdw6            = c6*VV;          
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            Vvdw12           = c12*VV;         
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;
            qq               = qH*jq;          
            r                = rsq21*rinv21;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            r                = rsq31*rinv31;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 12*n0;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}


