// For computing micro gemm for 8<M<16, N<8 and K, on the Posk-K computer.
// Parameger registers are pointers:
// x0: =>[M,N,K] TODO: is this long unsigned int?
// x1: =>[ALPHA,BETA]
// x2: =>A
// x3: =>LDA. Packed?
// x4: =>B
// x5: =>LDB. Packed?
// x6: =>C
// x7: =>LDC
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_dadddott16x8scalable
	.p2align	2
_dadddott16x8scalable:	// double-add-dot-16x8-scalable
	.cfi_startproc
	ldp	x9, x10, [x0], #16	// shape M and N.
	ldr	x8, [x0]	// loop parameter K.
	sub	x9, x9, #8	// X9 = M-8 TODO: determine by sign and branch to 8x8?
	orr	x15, xzr, #8	// use full space in first half of column dimension.
// Here only SVE+512bit is supported.
// Should determine Z-length to extend to other architectures. 
	whilelo p0.d, xzr, x15	// predicator in column first half
	whilelo	p1.d, xzr, x9	// ... in column last half
//	whilelo	p2.d, xzr, x10	// ... in row
// Zeroize C registers.
	mov	z13.d, xzr	// TODO: is this really equivalent to fp zero?
	mov	z14.d, xzr
	mov	z15.d, xzr
	mov	z16.d, xzr
	mov	z17.d, xzr
	mov	z18.d, xzr
	mov	z19.d, xzr
	mov	z20.d, xzr
	mov	z21.d, xzr
	mov	z22.d, xzr
	mov	z23.d, xzr
	mov	z24.d, xzr
	mov	z25.d, xzr
	mov	z26.d, xzr
	mov	z27.d, xzr
	mov	z28.d, xzr
	cmp	x8, #0
	cbz	WRITE_MEM	// branch directly to store.
K_START:
// Load column from A.
	ld1w	z0.d, p0/z, [x2]
	ld1w	z1.d, p1/z, [x2, #64]
	madd	x2, x2, x3, #8	// move forward
// Apply B columns.
	ldr	d0, [x4], #8	// col 1 of B, row 1 of C
	mov	z3.d, d0
	madd	z13.d, z13.d, z0.d, z3.d
	madd	z14.d, z14.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 2
	mov	z3.d, d0
	madd	z15.d, z15.d, z0.d, z3.d
	madd	z16.d, z16.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 3
	mov	z3.d, d0
	madd	z17.d, z17.d, z0.d, z3.d
	madd	z18.d, z18.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 4
	mov	z3.d, d0
	madd	z19.d, z19.d, z0.d, z3.d
	madd	z20.d, z20.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 5
	mov	z3.d, d0
	madd	z21.d, z21.d, z0.d, z3.d
	madd	z22.d, z22.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 6
	mov	z3.d, d0
	madd	z23.d, z23.d, z0.d, z3.d
	madd	z24.d, z24.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 7
	mov	z3.d, d0
	madd	z25.d, z25.d, z0.d, z3.d
	madd	z26.d, z26.d, z1.d, z3.d
	ldr	d0, [x4], #8	// 8
	mov	z3.d, d0
	madd	z27.d, z27.d, z0.d, z3.d
	madd	z28.d, z28.d, z1.d, z3.d
// TODO: optimize for packed situations.
	sub	x4, x4, #64
	madd	x4, x4, x5, #8
	subs	x8, x8, #1
	b.ne	K_START	// Next column of A & B.
WRITE_MEM:
// Z1&Z2 are now constant coefficients.
	ldr	d0, [x1], #8	// alpha
	mov	z0.d, d0
	ldr	d1, [x1]	// beta
	mov	z1.d, d1
// (R)Write data back to C memory.
	ld1w	z3.d, p0/z, [x6]	// 1
	mul	z3.d, z3.d, z1.d	// multiply by beta
	madd	z3.d, z13.d, z0.d	// colC = colC + alpha*DcolC
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z14.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8	// move pointer LDC doubles right
	ld1w	z3.d, p0/z, [x6]	// 2
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z15.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z16.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 3
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z17.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z18.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 4
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z19.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z20.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 5
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z21.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z22.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 6
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z23.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z24.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 7
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z25.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z26.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8
	ld1w	z3.d, p0/z, [x6]	// 8
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z27.d, z0.d
	st1w	z3.d, p0/z, [x6]
	ld1w	z3.d, p1/z, [x6, #64]
	mul	z3.d, z3.d, z1.d
	madd	z3.d, z28.d, z0.d
	st1w	z3.d, p1/z, [x6, #64]
	madd	x6, x6, x7, #8	// final increment is not used.
// PRELOAD_NEXT:
// TODO: do preloading here.
// End of calculation.
SUB_END:
	ret
	.cfi_endproc
// End of function.
.subsection_via_symbols.

