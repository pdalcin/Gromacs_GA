;#
;# $Id$
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as



.globl nb_kernel313_x86_64_sse
.globl _nb_kernel313_x86_64_sse
nb_kernel313_x86_64_sse:	
_nb_kernel313_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb313_fshift,           16
.equiv          nb313_gid,              24
.equiv          nb313_pos,              32
.equiv          nb313_faction,          40
.equiv          nb313_charge,           48
.equiv          nb313_p_facel,          56
.equiv          nb313_argkrf,           64
.equiv          nb313_argcrf,           72
.equiv          nb313_Vc,               80
.equiv          nb313_type,             88
.equiv          nb313_p_ntype,          96
.equiv          nb313_vdwparam,         104
.equiv          nb313_Vvdw,             112
.equiv          nb313_p_tabscale,       120
.equiv          nb313_VFtab,            128
.equiv          nb313_invsqrta,         136
.equiv          nb313_dvda,             144
.equiv          nb313_p_gbtabscale,     152
.equiv          nb313_GBtab,            160
.equiv          nb313_p_nthreads,       168
.equiv          nb313_count,            176
.equiv          nb313_mtx,              184
.equiv          nb313_outeriter,        192
.equiv          nb313_inneriter,        200
.equiv          nb313_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb313_ixO,              0
.equiv          nb313_iyO,              16
.equiv          nb313_izO,              32
.equiv          nb313_ixH1,             48
.equiv          nb313_iyH1,             64
.equiv          nb313_izH1,             80
.equiv          nb313_ixH2,             96
.equiv          nb313_iyH2,             112
.equiv          nb313_izH2,             128
.equiv          nb313_ixM,              144
.equiv          nb313_iyM,              160
.equiv          nb313_izM,              176
.equiv          nb313_iqM,              192
.equiv          nb313_iqH,              208
.equiv          nb313_epsH1,            224
.equiv          nb313_epsH2,            240
.equiv          nb313_epsM,             256
.equiv          nb313_dxH1,             272
.equiv          nb313_dyH1,             288
.equiv          nb313_dzH1,             304
.equiv          nb313_dxH2,             320
.equiv          nb313_dyH2,             336
.equiv          nb313_dzH2,             352
.equiv          nb313_dxM,              368
.equiv          nb313_dyM,              384
.equiv          nb313_dzM,              400
.equiv          nb313_qqM,              416
.equiv          nb313_qqH,              432
.equiv          nb313_rinvH1,           448
.equiv          nb313_rinvH2,           464
.equiv          nb313_rinvM,            480
.equiv          nb313_rH1,              496
.equiv          nb313_rH2,              512
.equiv          nb313_rM,               528
.equiv          nb313_tsc,              544
.equiv          nb313_two,              560
.equiv          nb313_c6,               576
.equiv          nb313_c12,              592
.equiv          nb313_six,              608
.equiv          nb313_twelve,           624
.equiv          nb313_vctot,            640
.equiv          nb313_Vvdwtot,          656
.equiv          nb313_fixO,             672
.equiv          nb313_fiyO,             688
.equiv          nb313_fizO,             704
.equiv          nb313_fixH1,            720
.equiv          nb313_fiyH1,            736
.equiv          nb313_fizH1,            752
.equiv          nb313_fixH2,            768
.equiv          nb313_fiyH2,            784
.equiv          nb313_fizH2,            800
.equiv          nb313_fixM,             816
.equiv          nb313_fiyM,             832
.equiv          nb313_fizM,             848
.equiv          nb313_fjx,              864
.equiv          nb313_fjy,              880
.equiv          nb313_fjz,              896
.equiv          nb313_half,             912
.equiv          nb313_three,            928
.equiv          nb313_is3,              944
.equiv          nb313_ii3,              948
.equiv          nb313_nri,              952
.equiv          nb313_iinr,             960
.equiv          nb313_jindex,           968
.equiv          nb313_jjnr,             976
.equiv          nb313_shift,            984
.equiv          nb313_shiftvec,         992
.equiv          nb313_facel,            1000
.equiv          nb313_innerjjnr,        1008
.equiv          nb313_ntia,             1016
.equiv          nb313_innerk,           1020
.equiv          nb313_n,                1024
.equiv          nb313_nn1,              1028
.equiv          nb313_nouter,           1032
.equiv          nb313_ninner,           1036
	push rbp
	mov  rbp, rsp
	push rbx

	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 1048		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb313_nouter], eax
	mov [rsp + nb313_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb313_nri], edi
	mov [rsp + nb313_iinr], rsi
	mov [rsp + nb313_jindex], rdx
	mov [rsp + nb313_jjnr], rcx
	mov [rsp + nb313_shift], r8
	mov [rsp + nb313_shiftvec], r9
	mov rsi, [rbp + nb313_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb313_facel], xmm0

	mov rax, [rbp + nb313_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb313_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb313_half], eax
	movss xmm1, [rsp + nb313_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + nb313_half],  xmm1
	movaps [rsp + nb313_two],  xmm2
	movaps [rsp + nb313_three],  xmm3
	movaps [rsp + nb313_six],  xmm4
	movaps [rsp + nb313_twelve],  xmm5
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb313_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb313_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb313_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb313_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb313_iqM], xmm3
	movaps [rsp + nb313_iqH], xmm4
	
	mov   rdx, [rbp + nb313_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb313_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb313_ntia], ecx		
.nb313_threadloop:
        mov   rsi, [rbp + nb313_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb313_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb313_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb313_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb313_n], eax
        mov [rsp + nb313_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb313_outerstart
        jmp .nb313_end
	
.nb313_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb313_nouter]
	mov [rsp + nb313_nouter], ebx

.nb313_outer:
	mov   rax, [rsp + nb313_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb313_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb313_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb313_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb313_pos]	;# rax = base of pos[]  
	mov   [rsp + nb313_ii3], ebx

	addss xmm3, [rax + rbx*4]  	;# ox
	addss xmm4, [rax + rbx*4 + 4]  ;# oy
	addss xmm5, [rax + rbx*4 + 8]  ;# oz
	addss xmm6, [rax + rbx*4 + 12] ;# h1x
	addss xmm7, [rax + rbx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [rsp + nb313_ixO], xmm3
	movaps [rsp + nb313_iyO], xmm4
	movaps [rsp + nb313_izO], xmm5
	movaps [rsp + nb313_ixH1], xmm6
	movaps [rsp + nb313_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [rax + rbx*4 + 20] ;# h1z
	addss xmm0, [rax + rbx*4 + 24] ;# h2x
	addss xmm1, [rax + rbx*4 + 28] ;# h2y
	addss xmm2, [rax + rbx*4 + 32] ;# h2z
	addss xmm3, [rax + rbx*4 + 36] ;# mx
	addss xmm4, [rax + rbx*4 + 40] ;# my
	addss xmm5, [rax + rbx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb313_izH1], xmm6
	movaps [rsp + nb313_ixH2], xmm0
	movaps [rsp + nb313_iyH2], xmm1
	movaps [rsp + nb313_izH2], xmm2
	movaps [rsp + nb313_ixM], xmm3
	movaps [rsp + nb313_iyM], xmm4
	movaps [rsp + nb313_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb313_vctot], xmm4
	movaps [rsp + nb313_Vvdwtot], xmm4
	movaps [rsp + nb313_fixO], xmm4
	movaps [rsp + nb313_fiyO], xmm4
	movaps [rsp + nb313_fizO], xmm4
	movaps [rsp + nb313_fixH1], xmm4
	movaps [rsp + nb313_fiyH1], xmm4
	movaps [rsp + nb313_fizH1], xmm4
	movaps [rsp + nb313_fixH2], xmm4
	movaps [rsp + nb313_fiyH2], xmm4
	movaps [rsp + nb313_fizH2], xmm4
	movaps [rsp + nb313_fixM], xmm4
	movaps [rsp + nb313_fiyM], xmm4
	movaps [rsp + nb313_fizM], xmm4
	
	mov   rax, [rsp + nb313_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb313_pos]
	mov   rdi, [rbp + nb313_faction]	
	mov   rax, [rsp + nb313_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb313_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb313_ninner]
	mov   [rsp + nb313_ninner], ecx
	add   edx, 0
	mov   [rsp + nb313_innerk], edx	;# number of innerloop atoms 
	jge   .nb313_unroll_loop
	jmp   .nb313_odd_inner
.nb313_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb313_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb313_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb313_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb313_iqM]
	mulps  xmm4, [rsp + nb313_iqH]

	movaps  [rsp + nb313_qqM], xmm3
	movaps  [rsp + nb313_qqH], xmm4
	
	mov rsi, [rbp + nb313_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb313_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb313_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	movlps xmm6, [rsi + r8*4]
	movlps xmm7, [rsi + r10*4]
	movhps xmm6, [rsi + r9*4]
	movhps xmm7, [rsi + r11*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movaps [rsp + nb313_c6], xmm4
	movaps [rsp + nb313_c12], xmm6

	mov rsi, [rbp + nb313_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
    
    ;# O interaction
    ;# copy to xmm3-xmm5
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    
    subps xmm3, [rsp + nb313_ixO]
    subps xmm4, [rsp + nb313_iyO]
    subps xmm5, [rsp + nb313_izO]
    
    movaps xmm13, xmm3
    movaps xmm14, xmm4
    movaps xmm15, xmm5
    
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm3, xmm4
	addps  xmm3, xmm5

    ;# calc 1/rsq
    rcpps xmm5, xmm3
    movaps xmm4, [rsp + nb313_two]
    mulps xmm3, xmm5
    subps xmm4, xmm3
    mulps xmm4, xmm5        ;# xmm4=rinvsq

    movaps xmm3, xmm4       ;# rinvsq
    mulps  xmm4, xmm4       ;# rinv4
    mulps  xmm4, xmm3       ;# rinv6
    movaps xmm5, xmm4     
    mulps  xmm5, xmm5       ;# rinv12
    mulps  xmm4, [rsp + nb313_c6]
    mulps  xmm5, [rsp + nb313_c12]
    movaps xmm6, xmm5
    subps  xmm6, xmm4  ;# Vvdw=vvdw12-vvdw6
    mulps  xmm4, [rsp + nb313_six]
    mulps  xmm5, [rsp + nb313_twelve]
    subps  xmm5, xmm4
    mulps  xmm3, xmm5   ;# fscal
    
    addps  xmm6, [rsp + nb313_Vvdwtot]
    movaps [rsp + nb313_Vvdwtot], xmm6
    
    mulps  xmm13, xmm3 ;# fx
    mulps  xmm14, xmm3 ;# fy
    mulps  xmm15, xmm3 ;# fz

    ;# save j force temporarily
    movaps [rsp + nb313_fjx], xmm13
    movaps [rsp + nb313_fjy], xmm14
    movaps [rsp + nb313_fjz], xmm15
    
    ;# increment i O force
    addps xmm13, [rsp + nb313_fixO]
    addps xmm14, [rsp + nb313_fiyO]
    addps xmm15, [rsp + nb313_fizO]
    movaps [rsp + nb313_fixO], xmm13
    movaps [rsp + nb313_fiyO], xmm14
    movaps [rsp + nb313_fizO], xmm15    
    ;# finished O LJ interaction.


    ;# do H1, H2, and M interactions in parallel.
    ;# xmm0-xmm2 still contain j coordinates.        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb313_ixH1]
    subps xmm1, [rsp + nb313_iyH1]
    subps xmm2, [rsp + nb313_izH1]
    subps xmm3, [rsp + nb313_ixH2]
    subps xmm4, [rsp + nb313_iyH2]
    subps xmm5, [rsp + nb313_izH2]
    subps xmm6, [rsp + nb313_ixM]
    subps xmm7, [rsp + nb313_iyM]
    subps xmm8, [rsp + nb313_izM]    

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps [rsp + nb313_dxH1], xmm0
	movaps [rsp + nb313_dyH1], xmm1
	movaps [rsp + nb313_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb313_dxH2], xmm3
	movaps [rsp + nb313_dyH2], xmm4
	movaps [rsp + nb313_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb313_dxM], xmm6
	movaps [rsp + nb313_dyM], xmm7
	movaps [rsp + nb313_dzM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for j atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb313_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb313_half]
	mulps   xmm9, xmm4  ;# rinvH1 
	mulps   xmm10, xmm4 ;# rinvH2
    mulps   xmm11, xmm4 ;# rinvM

	movaps  [rsp + nb313_rinvH1], xmm9
	movaps  [rsp + nb313_rinvH2], xmm10
	movaps  [rsp + nb313_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb313_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
        
    mov  rsi, [rbp + nb313_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb313_epsH1], xmm0
    movaps    [rsp + nb313_epsH2], xmm3
    movaps    [rsp + nb313_epsM], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb313_epsH1]
    movaps xmm13, [rsp + nb313_epsH2]
    movaps xmm14, [rsp + nb313_epsM]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb313_qqH]
    movaps xmm13, [rsp + nb313_qqM]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm12
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm12
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb313_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb313_vctot], xmm1
    
    movaps xmm10, [rsp + nb313_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
    
    movd eax, mm0
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

	;# move j forces to local temp variables 
    mov rdi, [rbp + nb313_faction]
    movlps xmm11, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm12, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm11, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm12, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm13, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm14, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm2,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm13, xmm2 ;# jza  -  jzb  -
    movlhps xmm14, xmm5 ;# jzc  -  jzd -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm11: jxa jya jxb jyb 
    ;# xmm12: jxc jyc jxd jyd
    ;# xmm13: jza jzb jzc jzd

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb313_rinvH1]
    mulps  xmm7, [rsp + nb313_rinvH2]
    mulps  xmm10, [rsp + nb313_rinvM]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb313_dxH1]
	mulps xmm1, [rsp + nb313_dyH1]
	mulps xmm2, [rsp + nb313_dzH1]
	mulps xmm3, [rsp + nb313_dxH2]
	mulps xmm4, [rsp + nb313_dyH2]
	mulps xmm5, [rsp + nb313_dzH2]
	mulps xmm6, [rsp + nb313_dxM]
	mulps xmm7, [rsp + nb313_dyM]
	mulps xmm8, [rsp + nb313_dzM]

    ;# fetch forces from O interaction
    movaps xmm14, [rsp + nb313_fjx]
    movaps xmm15, [rsp + nb313_fjy]
    addps  xmm13, [rsp + nb313_fjz]

    addps xmm14, xmm0
    addps xmm15, xmm1
    addps xmm13,  xmm2
    addps xmm0, [rsp + nb313_fixH1]
    addps xmm1, [rsp + nb313_fiyH1]
    addps xmm2, [rsp + nb313_fizH1]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb313_fixH2]
    addps xmm4, [rsp + nb313_fiyH2]
    addps xmm5, [rsp + nb313_fizH2]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb313_fixM]
    addps xmm7, [rsp + nb313_fiyM]
    addps xmm8, [rsp + nb313_fizM]

    movaps [rsp + nb313_fixH1], xmm0
    movaps [rsp + nb313_fiyH1], xmm1
    movaps [rsp + nb313_fizH1], xmm2
    movaps [rsp + nb313_fixH2], xmm3
    movaps [rsp + nb313_fiyH2], xmm4
    movaps [rsp + nb313_fizH2], xmm5
    movaps [rsp + nb313_fixM], xmm6
    movaps [rsp + nb313_fiyM], xmm7
    movaps [rsp + nb313_fizM], xmm8
    
    ;# xmm14 = fjx
    ;# xmm15 = fjy
    ;# xmm13 = fjz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm11
    movhps [rdi + rbx*4], xmm11
    movlps [rdi + rcx*4], xmm12
    movhps [rdi + rdx*4], xmm12
    movss  [rdi + rax*4 + 8], xmm13
    movss  [rdi + rcx*4 + 8], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 8], xmm13
    movss  [rdi + rdx*4 + 8], xmm14
    
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb313_innerk],  4
	jl    .nb313_odd_inner
	jmp   .nb313_unroll_loop
.nb313_odd_inner:	
	add dword ptr [rsp + nb313_innerk],  4
	jnz   .nb313_odd_loop
	jmp   .nb313_updateouterdata
.nb313_odd_loop:
	mov   rdx, [rsp + nb313_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb313_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb313_iqM]
	mov rsi, [rbp + nb313_charge] 
	movhps xmm4, [rsp + nb313_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb313_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb313_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb313_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb313_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb313_c6], xmm6
	movaps [rsp + nb313_c12], xmm7
	
	mov rsi, [rbp + nb313_pos]
	lea rax, [rax + rax*2]  

	movss xmm0, [rsp + nb313_ixO]
	movss xmm1, [rsp + nb313_iyO]
	movss xmm2, [rsp + nb313_izO]
	movss xmm3, [rsp + nb313_ixH1]
	movss xmm4, [rsp + nb313_iyH1]
	movss xmm5, [rsp + nb313_izH1]
	unpcklps xmm0, [rsp + nb313_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm1, [rsp + nb313_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm2, [rsp + nb313_izH2]	;# izO izH2 - -
	unpcklps xmm3, [rsp + nb313_ixM] 	;# ixH1 ixM - -
	unpcklps xmm4, [rsp + nb313_iyM]  	;# iyH1 iyM - -
	unpcklps xmm5, [rsp + nb313_izM]	;# izH1 izM - -
	unpcklps xmm0, xmm3  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm1, xmm4 	;# same for y
	unpcklps xmm2, xmm5 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# use M distances for storage
	movaps [rsp + nb313_dxM], xmm3
	movaps [rsp + nb313_dyM], xmm4
	movaps [rsp + nb313_dzM], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb313_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb313_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv
	
	movaps [rsp + nb313_rinvM], xmm0
	mulps  xmm4, xmm0  	;# r
	 
	mulps xmm4, [rsp + nb313_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm7, mm7
	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	mov  rsi, [rbp + nb313_VFtab]
    	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7

	xorps  xmm5, xmm5
	movlps xmm3, [rsi + rcx*4]	;# data: Y3 F3  -  - 
	movhps xmm5, [rsi + rbx*4]	;# data:  0  0 Y2 F2
	movhps xmm3, [rsi + rdx*4]      ;# data: Y3 F3 Y4 F4 

	movaps xmm4, xmm5		;# data:  0  0 Y2 F2 
	shufps xmm4, xmm3, 0x88		;# data:  0 Y2 Y3 Y3
	shufps xmm5, xmm3, 0xDD	        ;# data:  0 F2 F3 F4 

	xorps  xmm7, xmm7
	movlps xmm3, [rsi + rcx*4 + 8]	;# data: G3 H3  -  - 
	movhps xmm7, [rsi + rbx*4 + 8]	;# data:  0  0 G2 H2
	movhps xmm3, [rsi + rdx*4 + 8]  ;# data: G3 H3 G4 H4 

	movaps xmm6, xmm7		;# data:  0  0 G2 H2 
	shufps xmm6, xmm3, 0x88		;# data:  0 G2 G3 G3
	shufps xmm7, xmm3, 0xDD	        ;# data:  0 H2 H3 H4 

	;# xmm4 =  0  Y2 Y3 Y4
	;# xmm5 =  0  F2 F3 F4
	;# xmm6 =  0  G2 G3 G4
	;# xmm7 =  0  H2 H3 H4
	
	;# coulomb table ready, in xmm4-xmm7      
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	mulps  xmm7, [rsp + nb313_two]   	;# two*Heps2 
	movaps xmm0, [rsp + nb313_qqM]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	mulps  xmm0, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and xmm0 fijC 
	;# increment vcoul - then we can get rid of mm5 
	addps  xmm5, [rsp + nb313_vctot]
	movaps [rsp + nb313_vctot], xmm5
	
	;# do nontable L-J  in first element only.
	movaps xmm2, [rsp + nb313_rinvM]
	mulss  xmm2, xmm2
	movaps xmm1, xmm2
	mulss  xmm1, xmm1
	mulss  xmm1, xmm2	;# xmm1=rinvsix
	xorps  xmm4, xmm4
	movss  xmm4, xmm1
	mulss  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulss  xmm1, [rsp + nb313_c6]
	mulss  xmm4, [rsp + nb313_c12]
	movaps xmm3, xmm4
	subss  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	mulss  xmm1, [rsp + nb313_six]
	mulss  xmm4, [rsp + nb313_twelve]
	subss  xmm4, xmm1
	addss  xmm3, [rsp + nb313_Vvdwtot]
	mulss  xmm4, [rsp + nb313_rinvM]
	;# add back coul stuff from memory, and work on all elements again
	mulps  xmm0, [rsp + nb313_tsc]
	subps  xmm4, xmm0
	movss [rsp + nb313_Vvdwtot], xmm3
	mulps  xmm4, [rsp + nb313_rinvM]	
	
	movaps xmm0, [rsp + nb313_dxM]
	movaps xmm1, [rsp + nb313_dyM]
	movaps xmm2, [rsp + nb313_dzM]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force)
	
	movss  xmm3, [rsp + nb313_fixO]	
	movss  xmm4, [rsp + nb313_fiyO]	
	movss  xmm5, [rsp + nb313_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb313_fixO], xmm3	
	movss  [rsp + nb313_fiyO], xmm4	
	movss  [rsp + nb313_fizO], xmm5	;# updated the O force now do the H's

	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2      
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb313_fixH1]
	addss  xmm4, [rsp + nb313_fiyH1]
	addss  xmm5, [rsp + nb313_fizH1]
	movss  [rsp + nb313_fixH1], xmm3	
	movss  [rsp + nb313_fiyH1], xmm4	
	movss  [rsp + nb313_fizH1], xmm5	;# updated the H1 force 

	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb313_fixH2]
	addss  xmm4, [rsp + nb313_fiyH2]
	addss  xmm5, [rsp + nb313_fizH2]
	movss  [rsp + nb313_fixH2], xmm3	
	movss  [rsp + nb313_fiyH2], xmm4	
	movss  [rsp + nb313_fizH2], xmm5	;# updated the H2 force 

	mov rdi, [rbp + nb313_faction]
	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [rsp + nb313_fixM]
	addss  xmm4, [rsp + nb313_fiyM]
	addss  xmm5, [rsp + nb313_fizM]
	movss  [rsp + nb313_fixM], xmm3	
	movss  [rsp + nb313_fiyM], xmm4	
	movss  [rsp + nb313_fizM], xmm5	;# updated the M force 

	;# the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
	
	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  xmm0, xmm3
	movaps  xmm1, xmm4
	movaps  xmm2, xmm5
	
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5
	unpcklps xmm0, xmm1 	;# x,y sum in xmm0, z sum in xmm2
	
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb313_innerk]
	jz    .nb313_updateouterdata
	jmp   .nb313_odd_loop
.nb313_updateouterdata:
	mov   ecx, [rsp + nb313_ii3]
	mov   rdi, [rbp + nb313_faction]
	mov   rsi, [rbp + nb313_fshift]
	mov   edx, [rsp + nb313_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb313_fixO]
	movaps xmm1, [rsp + nb313_fiyO]
	movaps xmm2, [rsp + nb313_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb313_fixH1]
	movaps xmm1, [rsp + nb313_fiyH1]
	movaps xmm2, [rsp + nb313_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb313_fixH2]
	movaps xmm1, [rsp + nb313_fiyH2]
	movaps xmm2, [rsp + nb313_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb313_fixM]
	movaps xmm1, [rsp + nb313_fiyM]
	movaps xmm2, [rsp + nb313_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 36]
	movss  xmm4, [rdi + rcx*4 + 40]
	movss  xmm5, [rdi + rcx*4 + 44]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 36], xmm3
	movss  [rdi + rcx*4 + 40], xmm4
	movss  [rdi + rcx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb313_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb313_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb313_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb313_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb313_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb313_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb313_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb313_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb313_n], esi
        jmp .nb313_outer
.nb313_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb313_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb313_end
        ;# non-zero, do one more workunit
        jmp   .nb313_threadloop
.nb313_end:
	mov eax, [rsp + nb313_nouter]
	mov ebx, [rsp + nb313_ninner]
	mov rcx, [rbp + nb313_outeriter]
	mov rdx, [rbp + nb313_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1048
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret

	
	



.globl nb_kernel313nf_x86_64_sse
.globl _nb_kernel313nf_x86_64_sse
nb_kernel313nf_x86_64_sse:	
_nb_kernel313nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb313nf_fshift,         16
.equiv          nb313nf_gid,            24
.equiv          nb313nf_pos,            32
.equiv          nb313nf_faction,        40
.equiv          nb313nf_charge,         48
.equiv          nb313nf_p_facel,        56
.equiv          nb313nf_argkrf,         64
.equiv          nb313nf_argcrf,         72
.equiv          nb313nf_Vc,             80
.equiv          nb313nf_type,           88
.equiv          nb313nf_p_ntype,        96
.equiv          nb313nf_vdwparam,       104
.equiv          nb313nf_Vvdw,           112
.equiv          nb313nf_p_tabscale,     120
.equiv          nb313nf_VFtab,          128
.equiv          nb313nf_invsqrta,       136
.equiv          nb313nf_dvda,           144
.equiv          nb313nf_p_gbtabscale,   152
.equiv          nb313nf_GBtab,          160
.equiv          nb313nf_p_nthreads,     168
.equiv          nb313nf_count,          176
.equiv          nb313nf_mtx,            184
.equiv          nb313nf_outeriter,      192
.equiv          nb313nf_inneriter,      200
.equiv          nb313nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb313nf_ixO,            0
.equiv          nb313nf_iyO,            16
.equiv          nb313nf_izO,            32
.equiv          nb313nf_ixH1,           48
.equiv          nb313nf_iyH1,           64
.equiv          nb313nf_izH1,           80
.equiv          nb313nf_ixH2,           96
.equiv          nb313nf_iyH2,           112
.equiv          nb313nf_izH2,           128
.equiv          nb313nf_ixM,            144
.equiv          nb313nf_iyM,            160
.equiv          nb313nf_izM,            176
.equiv          nb313nf_iqM,            192
.equiv          nb313nf_iqH,            208
.equiv          nb313nf_qqM,            224
.equiv          nb313nf_qqH,            240
.equiv          nb313nf_rinvH1,         256
.equiv          nb313nf_rinvH2,         272
.equiv          nb313nf_rinvM,          288
.equiv          nb313nf_rH1,            304
.equiv          nb313nf_rH2,            320
.equiv          nb313nf_rM,             336
.equiv          nb313nf_tsc,            352
.equiv          nb313nf_two,            368
.equiv          nb313nf_c6,             384
.equiv          nb313nf_c12,            400
.equiv          nb313nf_vctot,          416
.equiv          nb313nf_Vvdwtot,        432
.equiv          nb313nf_half,           448
.equiv          nb313nf_three,          464
.equiv          nb313nf_is3,            480
.equiv          nb313nf_ii3,            484
.equiv          nb313nf_nri,            488
.equiv          nb313nf_iinr,           496
.equiv          nb313nf_jindex,         504
.equiv          nb313nf_jjnr,           512
.equiv          nb313nf_shift,          520
.equiv          nb313nf_shiftvec,       528
.equiv          nb313nf_facel,          536
.equiv          nb313nf_innerjjnr,      544
.equiv          nb313nf_ntia,           552
.equiv          nb313nf_innerk,         556
.equiv          nb313nf_n,              560
.equiv          nb313nf_nn1,            564
.equiv          nb313nf_nouter,         568
.equiv          nb313nf_ninner,         572
	push rbp
	mov  rbp, rsp
	push rbx

	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 584		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb313nf_nouter], eax
	mov [rsp + nb313nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb313nf_nri], edi
	mov [rsp + nb313nf_iinr], rsi
	mov [rsp + nb313nf_jindex], rdx
	mov [rsp + nb313nf_jjnr], rcx
	mov [rsp + nb313nf_shift], r8
	mov [rsp + nb313nf_shiftvec], r9
	mov rsi, [rbp + nb313nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb313nf_facel], xmm0

	mov rax, [rbp + nb313nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb313nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb313nf_half], eax
	movss xmm1, [rsp + nb313nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb313nf_half],  xmm1
	movaps [rsp + nb313nf_two],  xmm2
	movaps [rsp + nb313nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb313nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb313nf_charge]
	movss xmm4, [rdx + rbx*4 + 4]	
	movss xmm3, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb313nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb313nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb313nf_iqM], xmm3
	movaps [rsp + nb313nf_iqH], xmm4
	
	mov   rdx, [rbp + nb313nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb313nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb313nf_ntia], ecx		
.nb313nf_threadloop:
        mov   rsi, [rbp + nb313nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb313nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb313nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb313nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb313nf_n], eax
        mov [rsp + nb313nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb313nf_outerstart
        jmp .nb313nf_end
	
.nb313nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb313nf_nouter]
	mov [rsp + nb313nf_nouter], ebx

.nb313nf_outer:
	mov   rax, [rsp + nb313nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb313nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb313nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb313nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb313nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb313nf_ii3], ebx

	addss xmm3, [rax + rbx*4]  	;# ox
	addss xmm4, [rax + rbx*4 + 4]  ;# oy
	addss xmm5, [rax + rbx*4 + 8]  ;# oz
	addss xmm6, [rax + rbx*4 + 12] ;# h1x
	addss xmm7, [rax + rbx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [rsp + nb313nf_ixO], xmm3
	movaps [rsp + nb313nf_iyO], xmm4
	movaps [rsp + nb313nf_izO], xmm5
	movaps [rsp + nb313nf_ixH1], xmm6
	movaps [rsp + nb313nf_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [rax + rbx*4 + 20] ;# h1z
	addss xmm0, [rax + rbx*4 + 24] ;# h2x
	addss xmm1, [rax + rbx*4 + 28] ;# h2y
	addss xmm2, [rax + rbx*4 + 32] ;# h2z
	addss xmm3, [rax + rbx*4 + 36] ;# mx
	addss xmm4, [rax + rbx*4 + 40] ;# my
	addss xmm5, [rax + rbx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb313nf_izH1], xmm6
	movaps [rsp + nb313nf_ixH2], xmm0
	movaps [rsp + nb313nf_iyH2], xmm1
	movaps [rsp + nb313nf_izH2], xmm2
	movaps [rsp + nb313nf_ixM], xmm3
	movaps [rsp + nb313nf_iyM], xmm4
	movaps [rsp + nb313nf_izM], xmm5
	
	;# clear vctot
	xorps xmm4, xmm4
	movaps [rsp + nb313nf_vctot], xmm4
	movaps [rsp + nb313nf_Vvdwtot], xmm4

	mov   rax, [rsp + nb313nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb313nf_pos]
	mov   rdi, [rbp + nb313nf_faction]	
	mov   rax, [rsp + nb313nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb313nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb313nf_ninner]
	mov   [rsp + nb313nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb313nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb313nf_unroll_loop
	jmp   .nb313nf_odd_inner
.nb313nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb313nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb313nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb313nf_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb313nf_iqM]
	mulps  xmm4, [rsp + nb313nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb313nf_qqM], xmm3
	movaps  [rsp + nb313nf_qqH], xmm4
	
	mov rsi, [rbp + nb313nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb313nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb313nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb313nf_c6], xmm4
	movaps [rsp + nb313nf_c12], xmm6

	mov rsi, [rbp + nb313nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb313nf_ixO]
	movaps xmm5, [rsp + nb313nf_iyO]
	movaps xmm6, [rsp + nb313nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb313nf_ixH1]
	movaps xmm5, [rsp + nb313nf_iyH1]
	movaps xmm6, [rsp + nb313nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb313nf_ixH2]
	movaps xmm4, [rsp + nb313nf_iyH2]
	movaps xmm5, [rsp + nb313nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	
	;# move ixM-izM to xmm2-xmm4  
	movaps xmm3, [rsp + nb313nf_iyM]
	movaps xmm4, [rsp + nb313nf_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [rsp + nb313nf_ixM]
	subps  xmm2, xmm0	

	;# square it 
	mulps xmm2,xmm2
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	addps xmm4, xmm3
	addps xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb313nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb313nf_half]
	movaps  [rsp + nb313nf_rinvH1], xmm0	;# rinvH1 in xmm4 
	mulps   xmm6, xmm0
	movaps  [rsp + nb313nf_rH1], xmm6

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb313nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb313nf_half]
	movaps  [rsp + nb313nf_rinvH2], xmm0	;# rinvH2 in xmm4 
	mulps   xmm5, xmm0
	movaps  [rsp + nb313nf_rH2], xmm5

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [rsp + nb313nf_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [rsp + nb313nf_half]
	movaps  [rsp + nb313nf_rinvM], xmm0	;# rinvM in xmm5 
	mulps   xmm4, xmm0
	movaps  [rsp + nb313nf_rM], xmm4	
	
	;# Do the O LJ-only interaction directly.	
	rcpps   xmm2, xmm7
	movaps  xmm1, [rsp + nb313nf_two]
	mulps   xmm7, xmm2
	subps   xmm1, xmm7
	mulps   xmm2, xmm1 ;# rinvsq 
	movaps  xmm0, xmm2
	mulps   xmm0, xmm2  	;# r4
	mulps   xmm0, xmm2 	;# r6
	movaps  xmm1, xmm0
	mulps   xmm1, xmm1  	;# r12
	mulps   xmm0, [rsp + nb313nf_c6]
	mulps   xmm1, [rsp + nb313nf_c12]
	movaps  xmm3, xmm1
	subps   xmm3, xmm0  	;# Vvdw12-Vvdw6
	addps   xmm3, [rsp + nb313nf_Vvdwtot]
	movaps  [rsp + nb313nf_Vvdwtot], xmm3

	;# Do H1 interaction
	mov  rsi, [rbp + nb313nf_VFtab]
	
	movaps xmm7, [rsp + nb313nf_rH1]
	mulps   xmm7, [rsp + nb313nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx
	
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7      
        
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	movaps xmm0, [rsp + nb313nf_qqH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb313nf_vctot]
	movaps [rsp + nb313nf_vctot], xmm5 

	;# Done with H1, do H2 interactions 
	movaps xmm7, [rsp + nb313nf_rH2]
	mulps   xmm7, [rsp + nb313nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2
	
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# shuffle 10001000
	shufps xmm5, xmm7, 221  ;# shuffle 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# shuf 10001000
	shufps xmm7, xmm3, 221  ;# shuf 11011101
	;# coulomb table ready, in xmm4-xmm7      
        
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	movaps xmm0, [rsp + nb313nf_qqH]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb313nf_vctot]
	movaps [rsp + nb313nf_vctot], xmm5 

	;# Done with H2, do M interactions 
	movaps xmm7, [rsp + nb313nf_rM]
	mulps   xmm7, [rsp + nb313nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4	;# mm6/mm7 contain lu indices 
	
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm4, mm7
	movlhps xmm3, xmm4
	
	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2
	
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7      
        
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	movaps xmm0, [rsp + nb313nf_qqM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb313nf_vctot]
	movaps [rsp + nb313nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb313nf_innerk],  4
	jl    .nb313nf_odd_inner
	jmp   .nb313nf_unroll_loop
.nb313nf_odd_inner:	
	add dword ptr [rsp + nb313nf_innerk],  4
	jnz   .nb313nf_odd_loop
	jmp   .nb313nf_updateouterdata
.nb313nf_odd_loop:
	mov   rdx, [rsp + nb313nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb313nf_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [rsp + nb313nf_iqM]
	mov rsi, [rbp + nb313nf_charge] 
	movhps xmm4, [rsp + nb313nf_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb313nf_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov rsi, [rbp + nb313nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb313nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb313nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb313nf_c6], xmm6
	movaps [rsp + nb313nf_c12], xmm7
	
	mov rsi, [rbp + nb313nf_pos]
	lea rax, [rax + rax*2]  

	movss xmm3, [rsp + nb313nf_ixO]
	movss xmm4, [rsp + nb313nf_iyO]
	movss xmm5, [rsp + nb313nf_izO]
	movss xmm0, [rsp + nb313nf_ixH1]
	movss xmm1, [rsp + nb313nf_iyH1]
	movss xmm2, [rsp + nb313nf_izH1]
	unpcklps xmm3, [rsp + nb313nf_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [rsp + nb313nf_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [rsp + nb313nf_izH2]	;# izO izH2 - -
	unpcklps xmm0, [rsp + nb313nf_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [rsp + nb313nf_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [rsp + nb313nf_izM]	;# izH1 izM - -
	unpcklps xmm3, xmm0  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm4, xmm1 	;# same for y
	unpcklps xmm5, xmm2 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb313nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb313nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv
	
	movaps [rsp + nb313nf_rinvM], xmm0
	mulps  xmm4, xmm0  	;# r
	 
	mulps xmm4, [rsp + nb313nf_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm3, mm6
	cvtpi2ps xmm7, mm7
	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	mov  rsi, [rbp + nb313nf_VFtab]
    	psrlq mm6, 32
	movd ebx, mm6
	movd ecx, mm7
	psrlq mm7, 32
	movd edx, mm7

	xorps  xmm5, xmm5
	movlps xmm3, [rsi + rcx*4]	;# data: Y3 F3  -  - 
	movhps xmm5, [rsi + rbx*4]	;# data:  0  0 Y2 F2
	movhps xmm3, [rsi + rdx*4]      ;# data: Y3 F3 Y4 F4 

	movaps xmm4, xmm5		;# data:  0  0 Y2 F2 
	shufps xmm4, xmm3, 0x88		;# data:  0 Y2 Y3 Y3
	shufps xmm5, xmm3, 0xDD	        ;# data:  0 F2 F3 F4 

	xorps  xmm7, xmm7
	movlps xmm3, [rsi + rcx*4 + 8]	;# data: G3 H3  -  - 
	movhps xmm7, [rsi + rbx*4 + 8]	;# data:  0  0 G2 H2
	movhps xmm3, [rsi + rdx*4 + 8]  ;# data: G3 H3 G4 H4 

	movaps xmm6, xmm7		;# data:  0  0 G2 H2 
	shufps xmm6, xmm3, 0x88		;# data:  0 G2 G3 G3
	shufps xmm7, xmm3, 0xDD	        ;# data:  0 H2 H3 H4 

	;# xmm4 =  0  Y2 Y3 Y4
	;# xmm5 =  0  F2 F3 F4
	;# xmm6 =  0  G2 G3 G4
	;# xmm7 =  0  H2 H3 H4
	
	;# coulomb table ready, in xmm4-xmm7      
	mulps  xmm6, xmm1   	;# xmm6=Geps 
	mulps  xmm7, xmm2   	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7   	;# xmm5=Fp        
	movaps xmm0, [rsp + nb313nf_qqM]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul 
	addps  xmm5, [rsp + nb313nf_vctot]
	movaps [rsp + nb313nf_vctot], xmm5
	
	;# do nontable L-J  in first element only.
	movaps xmm2, [rsp + nb313nf_rinvM]
	mulss  xmm2, xmm2
	movaps xmm1, xmm2
	mulss  xmm1, xmm1
	mulss  xmm1, xmm2	;# xmm1=rinvsix
	xorps  xmm4, xmm4
	movss  xmm4, xmm1
	mulss  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulss  xmm1, [rsp + nb313nf_c6]
	mulss  xmm4, [rsp + nb313nf_c12]
	movaps xmm3, xmm4
	subss  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6
	addss  xmm3, [rsp + nb313nf_Vvdwtot]
	movss [rsp + nb313nf_Vvdwtot], xmm3

	dec dword ptr [rsp + nb313nf_innerk]
	jz    .nb313nf_updateouterdata
	jmp   .nb313nf_odd_loop
.nb313nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb313nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb313nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb313nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb313nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb313nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb313nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb313nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb313nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb313nf_n], esi
        jmp .nb313nf_outer
.nb313nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb313nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb313nf_end
        ;# non-zero, do one more workunit
        jmp   .nb313nf_threadloop
.nb313nf_end:
	mov eax, [rsp + nb313nf_nouter]
	mov ebx, [rsp + nb313nf_ninner]
	mov rcx, [rbp + nb313nf_outeriter]
	mov rdx, [rbp + nb313nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 584
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret


	
	
