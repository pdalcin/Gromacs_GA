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






.globl nb_kernel110_x86_64_sse2
.globl _nb_kernel110_x86_64_sse2
nb_kernel110_x86_64_sse2:	
_nb_kernel110_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb110_fshift,           16
.equiv          nb110_gid,              24
.equiv          nb110_pos,              32
.equiv          nb110_faction,          40
.equiv          nb110_charge,           48
.equiv          nb110_p_facel,          56
.equiv          nb110_argkrf,           64
.equiv          nb110_argcrf,           72
.equiv          nb110_Vc,               80
.equiv          nb110_type,             88
.equiv          nb110_p_ntype,          96
.equiv          nb110_vdwparam,         104
.equiv          nb110_Vvdw,             112
.equiv          nb110_p_tabscale,       120
.equiv          nb110_VFtab,            128
.equiv          nb110_invsqrta,         136
.equiv          nb110_dvda,             144
.equiv          nb110_p_gbtabscale,     152
.equiv          nb110_GBtab,            160
.equiv          nb110_p_nthreads,       168
.equiv          nb110_count,            176
.equiv          nb110_mtx,              184
.equiv          nb110_outeriter,        192
.equiv          nb110_inneriter,        200
.equiv          nb110_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb110_ix,               0
.equiv          nb110_iy,               16
.equiv          nb110_iz,               32
.equiv          nb110_iq,               48
.equiv          nb110_dx,               64
.equiv          nb110_dy,               80
.equiv          nb110_dz,               96
.equiv          nb110_c6,               112
.equiv          nb110_c12,              128
.equiv          nb110_six,              144
.equiv          nb110_twelve,           160
.equiv          nb110_vctot,            176
.equiv          nb110_Vvdwtot,          192
.equiv          nb110_fix,              208
.equiv          nb110_fiy,              224
.equiv          nb110_fiz,              240
.equiv          nb110_half,             256
.equiv          nb110_three,            272
.equiv          nb110_is3,              288
.equiv          nb110_ii3,              292
.equiv          nb110_nri,              296
.equiv          nb110_iinr,             304
.equiv          nb110_jindex,           312
.equiv          nb110_jjnr,             320
.equiv          nb110_shift,            328
.equiv          nb110_shiftvec,         336
.equiv          nb110_facel,            344
.equiv          nb110_innerjjnr,        352
.equiv          nb110_ntia,             360
.equiv          nb110_innerk,           364
.equiv          nb110_n,                368
.equiv          nb110_nn1,              372
.equiv          nb110_ntype,            376
.equiv          nb110_nouter,           380
.equiv          nb110_ninner,           384
	push rbp
	mov  rbp, rsp
	push rbx
	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 408		;# local variable stack space (n*16+8)
	
	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb110_nouter], eax
	mov [rsp + nb110_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb110_nri], edi
	mov [rsp + nb110_iinr], rsi
	mov [rsp + nb110_jindex], rdx
	mov [rsp + nb110_jjnr], rcx
	mov [rsp + nb110_shift], r8
	mov [rsp + nb110_shiftvec], r9
	mov rdi, [rbp + nb110_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb110_ntype], edi
	mov rsi, [rbp + nb110_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb110_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb110_half], eax
	mov [rsp + nb110_half+4], ebx
	movsd xmm1, [rsp + nb110_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb110_half], xmm1
	movapd [rsp + nb110_three], xmm3
	movapd [rsp + nb110_six], xmm4
	movapd [rsp + nb110_twelve], xmm5

.nb110_threadloop:
        mov   rsi, [rbp + nb110_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb110_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb110_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb110_n], eax
        mov [rsp + nb110_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110_outerstart
	jmp .nb110_end
	
.nb110_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb110_nouter]
	mov [rsp + nb110_nouter], ebx

.nb110_outer:
	mov   rax, [rsp + nb110_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb110_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb110_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb110_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb110_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb110_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb110_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb110_ntype]
    	shl   edx, 1
    	mov   [rsp + nb110_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb110_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb110_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb110_ix], xmm0
	movapd [rsp + nb110_iy], xmm1
	movapd [rsp + nb110_iz], xmm2

	mov   [rsp + nb110_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm12, xmm12
	movapd [rsp + nb110_Vvdwtot], xmm12
	movapd xmm13, xmm12
	movapd xmm14, xmm12
	movapd xmm15, xmm12
		
	mov   rax, [rsp + nb110_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb110_pos]
	mov   rdi, [rbp + nb110_faction]	
	mov   rax, [rsp + nb110_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb110_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb110_ninner]
	mov   [rsp + nb110_ninner], ecx
	add   edx, 0
	mov   [rsp + nb110_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110_unroll_loop
	jmp   .nb110_checksingle
.nb110_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	add qword ptr [rsp + nb110_innerjjnr],  8	;# advance pointer (unrolled 2) 
	
	mov rsi, [rbp + nb110_pos]       ;# base of pos[] 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb110_ix]
	subpd xmm5, [rsp + nb110_iy]
	subpd xmm6, [rsp + nb110_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	mov rsi, [rbp + nb110_charge]    ;# base of charge[] 

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mov rdi, [rbp + nb110_type]
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	movlpd xmm3, [rsi + r8*8]
	movhpd xmm3, [rsi + r9*8]

	cvtpd2ps xmm5, xmm4	
	mov r8d, [rdi + r8*4]
	mov r9d, [rdi + r9*4]
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	shl r8d, 1
	shl r9d, 1

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	mov edi, [rsp + nb110_ntia]
	add r8d, edi
	add r9d, edi
	movapd xmm1, [rsp + nb110_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	mov rdi, [rbp + nb110_vdwparam]
	movapd xmm0, [rsp + nb110_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 
	mulpd xmm3, [rsp + nb110_iq]		;# qq 

    movlpd xmm6, [rdi + r8*8]
    movlpd xmm7, [rdi + r8*8 + 8]
    
	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb110_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
    movhpd xmm6, [rdi + r9*8]
    movhpd xmm7, [rdi + r9*8 + 8]

	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm1, xmm6
	mulpd  xmm2, xmm7
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb110_Vvdwtot]
	mulpd  xmm1, [rsp + nb110_six]
	mulpd  xmm2, [rsp + nb110_twelve]
	subpd  xmm2, xmm1
	addpd  xmm2, xmm3
	mulpd  xmm4, xmm2	;# xmm4=total fscal 
	addpd  xmm12, xmm3  ;# add to vctot
	mov    rdi, [rbp + nb110_faction]

	;# the fj's - start by accumulating forces from memory 
	movlpd xmm6, [rdi + rax*8]
	movlpd xmm7, [rdi + rax*8 + 8]
	movlpd xmm8, [rdi + rax*8 + 16]
	movhpd xmm6, [rdi + rbx*8]
	movhpd xmm7, [rdi + rbx*8 + 8]
	movhpd xmm8, [rdi + rbx*8 + 16]

	movapd [rsp + nb110_Vvdwtot], xmm5

	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4
    
	addpd xmm6, xmm9
	addpd xmm7, xmm10
	addpd xmm8, xmm11

	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11

	movlpd [rdi + rax*8], xmm6
	movlpd [rdi + rax*8 + 8], xmm7
	movlpd [rdi + rax*8 + 16], xmm8
	movhpd [rdi + rbx*8], xmm6
	movhpd [rdi + rbx*8 + 8], xmm7
	movhpd [rdi + rbx*8 + 16], xmm8
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb110_innerk],  2
	jl    .nb110_checksingle
	jmp   .nb110_unroll_loop	
.nb110_checksingle:
	mov   edx, [rsp + nb110_innerk]
	and   edx, 1
	jnz    .nb110_dosingle
	jmp    .nb110_updateouterdata
.nb110_dosingle:
	mov rsi, [rbp + nb110_charge]
	mov rdi, [rbp + nb110_pos]
	mov rcx, [rsp + nb110_innerjjnr]
	mov   eax, [rcx]
	
	mov rsi, [rbp + nb110_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]
    mulsd xmm3, [rsp + nb110_iq] ;# qq

	mov rsi, [rbp + nb110_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb110_vdwparam]
	shl r8d, 1
	mov edi, [rsp + nb110_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]	    ;# c6
	movsd xmm6, [rsi + r8*8 + 8]	;# c12
	movapd [rsp + nb110_c6], xmm4
	movapd [rsp + nb110_c12], xmm6
	
	mov rsi, [rbp + nb110_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm4-xmm6
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]
	
	;# calc dr 
	subsd xmm4, [rsp + nb110_ix]
	subsd xmm5, [rsp + nb110_iy]
	subsd xmm6, [rsp + nb110_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb110_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb110_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm1, [rsp + nb110_c6]
	mulsd  xmm2, [rsp + nb110_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb110_Vvdwtot]
	mulsd  xmm1, [rsp + nb110_six]
	mulsd  xmm2, [rsp + nb110_twelve]
	subsd  xmm2, xmm1
	addsd  xmm2, xmm3
	mulsd  xmm4, xmm2	;# xmm4=total fscal 
	addsd  xmm12, xmm3   ;# add to vctot

	movsd [rsp + nb110_Vvdwtot], xmm5

	mov    rdi, [rbp + nb110_faction]
	mulsd  xmm9, xmm4
	mulsd  xmm10, xmm4
	mulsd  xmm11, xmm4
    
	;# now update f_i 
	addsd  xmm13, xmm9
	addsd  xmm14, xmm10
	addsd  xmm15, xmm11
	;# the fj's - start by accumulating forces from memory 
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
	
.nb110_updateouterdata:
	mov   ecx, [rsp + nb110_ii3]
	mov   rdi, [rbp + nb110_faction]
	mov   rsi, [rbp + nb110_fshift]
	mov   edx, [rsp + nb110_is3]

	;# accumulate i forces in xmm13, xmm14, xmm15
	movhlps xmm3, xmm13
	movhlps xmm4, xmm14
	movhlps xmm5, xmm15
	addsd  xmm13, xmm3
	addsd  xmm14, xmm4
	addsd  xmm15, xmm5 ;# sum is in low xmm13-xmm15

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [rsi + rdx*8]
	movsd  xmm4, [rsi + rdx*8 + 8]
	movsd  xmm5, [rsi + rdx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb110_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb110_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movhlps xmm6, xmm12
	addsd  xmm12, xmm6	;# low xmm12 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb110_Vc]
	addsd xmm12, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm12
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb110_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb110_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb110_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb110_n], esi
        jmp .nb110_outer
.nb110_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb110_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110_end
        ;# non-zero, do one more workunit
        jmp   .nb110_threadloop
.nb110_end:
	mov eax, [rsp + nb110_nouter]
	mov ebx, [rsp + nb110_ninner]
	mov rcx, [rbp + nb110_outeriter]
	mov rdx, [rbp + nb110_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 408
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret






.globl nb_kernel110nf_x86_64_sse2
.globl _nb_kernel110nf_x86_64_sse2
nb_kernel110nf_x86_64_sse2:	
_nb_kernel110nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb110nf_fshift,         16
.equiv          nb110nf_gid,            24
.equiv          nb110nf_pos,            32
.equiv          nb110nf_faction,        40
.equiv          nb110nf_charge,         48
.equiv          nb110nf_p_facel,        56
.equiv          nb110nf_argkrf,         64
.equiv          nb110nf_argcrf,         72
.equiv          nb110nf_Vc,             80
.equiv          nb110nf_type,           88
.equiv          nb110nf_p_ntype,        96
.equiv          nb110nf_vdwparam,       104
.equiv          nb110nf_Vvdw,           112
.equiv          nb110nf_p_tabscale,     120
.equiv          nb110nf_VFtab,          128
.equiv          nb110nf_invsqrta,       136
.equiv          nb110nf_dvda,           144
.equiv          nb110nf_p_gbtabscale,   152
.equiv          nb110nf_GBtab,          160
.equiv          nb110nf_p_nthreads,     168
.equiv          nb110nf_count,          176
.equiv          nb110nf_mtx,            184
.equiv          nb110nf_outeriter,      192
.equiv          nb110nf_inneriter,      200
.equiv          nb110nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb110nf_ix,             0
.equiv          nb110nf_iy,             16
.equiv          nb110nf_iz,             32
.equiv          nb110nf_iq,             48
.equiv          nb110nf_c6,             64
.equiv          nb110nf_c12,            80
.equiv          nb110nf_vctot,          96
.equiv          nb110nf_Vvdwtot,        112
.equiv          nb110nf_half,           128
.equiv          nb110nf_three,          144
.equiv          nb110nf_is3,            160
.equiv          nb110nf_ii3,            164
.equiv          nb110nf_nri,            168
.equiv          nb110nf_iinr,           176
.equiv          nb110nf_jindex,         184
.equiv          nb110nf_jjnr,           192
.equiv          nb110nf_shift,          200
.equiv          nb110nf_shiftvec,       208
.equiv          nb110nf_facel,          216
.equiv          nb110nf_innerjjnr,      224
.equiv          nb110nf_ntia,           232
.equiv          nb110nf_innerk,         236
.equiv          nb110nf_n,              240
.equiv          nb110nf_nn1,            244
.equiv          nb110nf_ntype,          248
.equiv          nb110nf_nouter,         252
.equiv          nb110nf_ninner,         256
	push rbp
	mov  rbp, rsp
	push rbx
	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 280		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb110nf_nouter], eax
	mov [rsp + nb110nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb110nf_nri], edi
	mov [rsp + nb110nf_iinr], rsi
	mov [rsp + nb110nf_jindex], rdx
	mov [rsp + nb110nf_jjnr], rcx
	mov [rsp + nb110nf_shift], r8
	mov [rsp + nb110nf_shiftvec], r9
	mov rdi, [rbp + nb110nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb110nf_ntype], edi
	mov rsi, [rbp + nb110nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb110nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb110nf_half], eax
	mov [rsp + nb110nf_half+4], ebx
	movsd xmm1, [rsp + nb110nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb110nf_half], xmm1
	movapd [rsp + nb110nf_three], xmm3

.nb110nf_threadloop:
        mov   rsi, [rbp + nb110nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb110nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb110nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb110nf_n], eax
        mov [rsp + nb110nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110nf_outerstart
        jmp .nb110nf_end

.nb110nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb110nf_nouter]
	mov [rsp + nb110nf_nouter], ebx

.nb110nf_outer:
	mov   rax, [rsp + nb110nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb110nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb110nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb110nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb110nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb110nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb110nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb110nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb110nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb110nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb110nf_ix], xmm0
	movapd [rsp + nb110nf_iy], xmm1
	movapd [rsp + nb110nf_iz], xmm2

	mov   [rsp + nb110nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb110nf_vctot], xmm4
	movapd [rsp + nb110nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb110nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb110nf_pos]
	mov   rax, [rsp + nb110nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb110nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb110nf_ninner]
	mov   [rsp + nb110nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb110nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110nf_unroll_loop
	jmp   .nb110nf_checksingle
.nb110nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb110nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb110nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb110nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm5, [rsp + nb110nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb110nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb110nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb110nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb110nf_c6], xmm4
	movapd [rsp + nb110nf_c12], xmm6
	
	mov rsi, [rbp + nb110nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb110nf_ix]
	movapd xmm5, [rsp + nb110nf_iy]
	movapd xmm6, [rsp + nb110nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb110nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb110nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm1, [rsp + nb110nf_c6]
	mulpd  xmm2, [rsp + nb110nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb110nf_Vvdwtot]
	addpd  xmm3, [rsp + nb110nf_vctot]
	movapd [rsp + nb110nf_vctot], xmm3
	movapd [rsp + nb110nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb110nf_innerk],  2
	jl    .nb110nf_checksingle
	jmp   .nb110nf_unroll_loop	
.nb110nf_checksingle:
	mov   edx, [rsp + nb110nf_innerk]
	and   edx, 1
	jnz   .nb110nf_dosingle
	jmp   .nb110nf_updateouterdata
.nb110nf_dosingle:
	mov rsi, [rbp + nb110nf_charge]
	mov rdi, [rbp + nb110nf_pos]
	mov rcx, [rsp + nb110nf_innerjjnr]
	mov   eax, [rcx]
	
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]

	movapd xmm5, [rsp + nb110nf_iq]
	mulsd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov rsi, [rbp + nb110nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb110nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb110nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	
	movapd [rsp + nb110nf_c6], xmm4
	movapd [rsp + nb110nf_c12], xmm6
	
	mov rsi, [rbp + nb110nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb110nf_ix]
	movapd xmm5, [rsp + nb110nf_iy]
	movapd xmm6, [rsp + nb110nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb110nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb110nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb110nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm1, [rsp + nb110nf_c6]
	mulsd  xmm2, [rsp + nb110nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb110nf_Vvdwtot]
	addsd  xmm3, [rsp + nb110nf_vctot]
	movlpd [rsp + nb110nf_vctot], xmm3
	movlpd [rsp + nb110nf_Vvdwtot], xmm5
	
.nb110nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb110nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb110nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb110nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb110nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb110nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb110nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb110nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb110nf_n], esi
        jmp .nb110nf_outer
.nb110nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb110nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110nf_end
        ;# non-zero, do one more workunit
        jmp   .nb110nf_threadloop
.nb110nf_end:
	mov eax, [rsp + nb110nf_nouter]
	mov ebx, [rsp + nb110nf_ninner]
	mov rcx, [rbp + nb110nf_outeriter]
	mov rdx, [rbp + nb110nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 280
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret	

