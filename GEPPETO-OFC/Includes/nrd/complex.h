#define ANSI  // should probably be somewhere else but I can't be bothered PB Jan 2004


#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#ifndef DCOMPLEX_DECLARE_T_
typedef struct DCOMPLEX {double r,i;} dcomplex;
#define _DCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

dcomplex Cadd(dcomplex a, dcomplex b);
dcomplex Csub(dcomplex a, dcomplex b);
dcomplex Cmul(dcomplex a, dcomplex b);
dcomplex Complex(float re, float im);
dcomplex Conjg(dcomplex z);
dcomplex Cdiv(dcomplex a, dcomplex b);
float Cabs(dcomplex z);
dcomplex Csqrt(dcomplex z);
dcomplex RCmul(float x, dcomplex a);

#else /* ANSI */
/* traditional - K&R */

dcomplex Cadd();
dcomplex Csub();
dcomplex Cmul();
dcomplex Complex();
dcomplex Conjg();
dcomplex Cdiv();
float Cabs();
dcomplex Csqrt();
dcomplex RCmul();

#endif /* ANSI */

#endif /* _NR_COMPLEX_H_ */
