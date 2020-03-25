// For computing micro gemm for W<M<2W, N<14 and K,
//   where W is 2*{SVE vector length in words}.
// Parameger registers are pointers:
// x0: =>[M,N,K]
// x1: =>[ALPHA,BETA]
// x2: =>A
// x3: LDA
// x4: =>B
// x5: LDB
// x6: =>C
// x7: LDC
	.arch	armv8.2-a+sve
	.text
	.global	_dmgemm_sve2wx4
	.type	_dmgemm_sve2wx4, %function
_dmgemm_sve2wx4:
	.cfi_startproc
	ldp	x9, x10, [x0], #16	// shape M and N.
	ldr	x8, [x0]	// loop parameter K.
	mov	x11, xzr
	mov	x12, #8	// double size
	incd	x11	// determine vector length, in doubles.
	ptrue	p0.d, all	// first half is all-true.
	whilelo	p1.d, x11, x9	// second half from M argument.
	fmov	d0, #1.0	// exact float 1.0.
	fmov	x14, d0	// do not conflict with sve registers.
// Register configuration:
// Z[30-31]: A columns
// Z29: B elements broadcasted
// Z28: multiplication buffer
// Z[0-27]: C change buffer
	dup	z0.d, #0
	dup	z1.d, #0
	dup	z2.d, #0
	dup	z3.d, #0
	dup	z4.d, #0
	dup	z5.d, #0
	dup	z6.d, #0
	dup	z7.d, #0
	dup	z8.d, #0
	dup	z9.d, #0
	dup	z10.d, #0
	dup	z11.d, #0
	dup	z12.d, #0
	dup	z13.d, #0
	dup	z14.d, #0
	dup	z15.d, #0
	dup	z16.d, #0
	dup	z17.d, #0
	dup	z18.d, #0
	dup	z19.d, #0
	dup	z20.d, #0
	dup	z21.d, #0
	dup	z22.d, #0
	dup	z23.d, #0
	dup	z24.d, #0
	dup	z25.d, #0
	dup	z26.d, #0
	dup	z27.d, #0
K_LOOP:
// Load columns from A.
	ld1d	z30.d, p0/z, [x2]
	ld1d	z31.d, p1/z, [x2, x11, lsl 3]	// second vector
	madd	x2, x3, x12, x2	// move forward
// Apply B columns.
	mov	x13, x10	// counter
	ld1rd	z29.d, p0/z, [x4]	// row L column 0
	fmul	z28.d, z30.d, z29.d
	fadd	z0.d, z0.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z1.d, z1.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #8]	// row L column 1
	fmul	z28.d, z30.d, z29.d
	fadd	z2.d, z2.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z3.d, z3.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #16]	// row L column 2
	fmul	z28.d, z30.d, z29.d
	fadd	z4.d, z4.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z5.d, z5.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #24]	// row L column 3
	fmul	z28.d, z30.d, z29.d
	fadd	z6.d, z6.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z7.d, z7.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #32]	// row L column 4
	fmul	z28.d, z30.d, z29.d
	fadd	z8.d, z8.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z9.d, z9.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #40]	// row L column 5
	fmul	z28.d, z30.d, z29.d
	fadd	z10.d, z10.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z11.d, z11.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #48]	// row L column 6
	fmul	z28.d, z30.d, z29.d
	fadd	z12.d, z12.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z13.d, z13.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #56]	// row L column 7
	fmul	z28.d, z30.d, z29.d
	fadd	z14.d, z14.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z15.d, z15.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #64]	// row L column 8
	fmul	z28.d, z30.d, z29.d
	fadd	z16.d, z16.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z17.d, z17.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #72]	// row L column 9
	fmul	z28.d, z30.d, z29.d
	fadd	z18.d, z18.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z19.d, z19.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #80]	// row L column 10
	fmul	z28.d, z30.d, z29.d
	fadd	z20.d, z20.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z21.d, z21.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #88]	// row L column 11
	fmul	z28.d, z30.d, z29.d
	fadd	z22.d, z22.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z23.d, z23.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #96]	// row L column 12
	fmul	z28.d, z30.d, z29.d
	fadd	z24.d, z24.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z25.d, z25.d, z28.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW
	ld1rd	z29.d, p0/z, [x4, #104]	// row L column 13
	fmul	z28.d, z30.d, z29.d
	fadd	z26.d, z26.d, z28.d
	fmul	z28.d, z31.d, z29.d
	fadd	z27.d, z27.d, z28.d
//	subs	x13, x13, #1
//	b	NEXT_ROW
NEXT_ROW:
	madd	x4, x5, x12, x4	// move forward
	subs	x8, x8, #1
	b.ne	K_LOOP	// next column / row.
WRITE_MEM:
// Override A and B buffers:
// z[30-31]: extended alpha and beta.
// z[28-29]: C memory buffer.
	ldr	x15, [x1]	// alpha, as 64-bits
	ld1rd	z30.d, p0/z, [x1]	// alpha, to the vector.
	ld1rd	z31.d, p0/z, [x1, #8]	// beta.
// (R&)Write data back to C memory.
	cmp	x14, x15
	b.eq	UNIT_ALPHA
// Non-unit alpha case.
// Scale all C change buffers.
	fmul	z0.d, z0.d, z30.d
	fmul	z1.d, z1.d, z30.d
	fmul	z2.d, z2.d, z30.d
	fmul	z3.d, z3.d, z30.d
	fmul	z4.d, z4.d, z30.d
	fmul	z5.d, z5.d, z30.d
	fmul	z6.d, z6.d, z30.d
	fmul	z7.d, z7.d, z30.d
	fmul	z8.d, z8.d, z30.d
	fmul	z9.d, z9.d, z30.d
	fmul	z10.d, z10.d, z30.d
	fmul	z11.d, z11.d, z30.d
	fmul	z12.d, z12.d, z30.d
	fmul	z13.d, z13.d, z30.d
	fmul	z14.d, z14.d, z30.d
	fmul	z15.d, z15.d, z30.d
	fmul	z16.d, z16.d, z30.d
	fmul	z17.d, z17.d, z30.d
	fmul	z18.d, z18.d, z30.d
	fmul	z19.d, z19.d, z30.d
	fmul	z20.d, z20.d, z30.d
	fmul	z21.d, z21.d, z30.d
	fmul	z22.d, z22.d, z30.d
	fmul	z23.d, z23.d, z30.d
	fmul	z24.d, z24.d, z30.d
	fmul	z25.d, z25.d, z30.d
	fmul	z26.d, z26.d, z30.d
	fmul	z27.d, z27.d, z30.d
// Unit alpha case.
UNIT_ALPHA:
//	mov	x10, x10	// x10 itself acts as counter.
	ld1d	z28.d, p0/z, [x6]	// column vector 0
	fmad	z28.d, p0/m, z31.d, z0.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z1.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 1
	fmad	z28.d, p0/m, z31.d, z2.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z3.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 2
	fmad	z28.d, p0/m, z31.d, z4.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z5.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 3
	fmad	z28.d, p0/m, z31.d, z6.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z7.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 4
	fmad	z28.d, p0/m, z31.d, z8.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z9.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 5
	fmad	z28.d, p0/m, z31.d, z10.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z11.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 6
	fmad	z28.d, p0/m, z31.d, z12.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z13.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 7
	fmad	z28.d, p0/m, z31.d, z14.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z15.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 8
	fmad	z28.d, p0/m, z31.d, z16.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z17.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 9
	fmad	z28.d, p0/m, z31.d, z18.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z19.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 10
	fmad	z28.d, p0/m, z31.d, z20.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z21.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 11
	fmad	z28.d, p0/m, z31.d, z22.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z23.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 12
	fmad	z28.d, p0/m, z31.d, z24.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z25.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM
	ld1d	z28.d, p0/z, [x6]	// column vector 13
	fmad	z28.d, p0/m, z31.d, z26.d
	st1d	z28.d, p0, [x6]
	ld1d	z29.d, p1/z, [x6, x11, lsl #3]
	fmad	z29.d, p1/m, z31.d, z27.d
	st1d	z29.d, p1, [x6, x11, lsl #3]
//	subs	x10, x10, #1
//	madd	x6, x7, x12, x6
//	b	END_WRITE_MEM
// End of computation.
END_WRITE_MEM:
	mov	x0, #0	// return normal.
	b	END_EXEC
END_ERROR:
	mov	x0, #1	// return error.
END_EXEC:
	ret
	.cfi_endproc
