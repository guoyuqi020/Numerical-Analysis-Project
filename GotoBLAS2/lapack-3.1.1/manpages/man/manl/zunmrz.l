.TH  ZUNMRZ 1 "February 2007" " LAPACK routine (version 3.1.1) " " LAPACK routine (version 3.1.1) " 
.SH NAME
ZUNMRZ - the general complex M-by-N matrix C with   SIDE = \(aqL\(aq SIDE = \(aqR\(aq TRANS = \(aqN\(aq
.SH SYNOPSIS
.TP 19
SUBROUTINE ZUNMRZ(
SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
WORK, LWORK, INFO )
.TP 19
.ti +4
CHARACTER
SIDE, TRANS
.TP 19
.ti +4
INTEGER
INFO, K, L, LDA, LDC, LWORK, M, N
.TP 19
.ti +4
COMPLEX*16
A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
.SH PURPOSE
ZUNMRZ overwrites the general complex M-by-N matrix C with
TRANS = \(aqC\(aq:      Q**H * C       C * Q**H
.br

where Q is a complex unitary matrix defined as the product of k
elementary reflectors
.br

      Q = H(1) H(2) . . . H(k)
.br

as returned by ZTZRZF. Q is of order M if SIDE = \(aqL\(aq and of order N
if SIDE = \(aqR\(aq.
.br

.SH ARGUMENTS
.TP 8
SIDE    (input) CHARACTER*1
= \(aqL\(aq: apply Q or Q**H from the Left;
.br
= \(aqR\(aq: apply Q or Q**H from the Right.
.TP 8
TRANS   (input) CHARACTER*1
.br
= \(aqN\(aq:  No transpose, apply Q;
.br
= \(aqC\(aq:  Conjugate transpose, apply Q**H.
.TP 8
M       (input) INTEGER
The number of rows of the matrix C. M >= 0.
.TP 8
N       (input) INTEGER
The number of columns of the matrix C. N >= 0.
.TP 8
K       (input) INTEGER
The number of elementary reflectors whose product defines
the matrix Q.
If SIDE = \(aqL\(aq, M >= K >= 0;
if SIDE = \(aqR\(aq, N >= K >= 0.
.TP 8
L       (input) INTEGER
The number of columns of the matrix A containing
the meaningful part of the Householder reflectors.
If SIDE = \(aqL\(aq, M >= L >= 0, if SIDE = \(aqR\(aq, N >= L >= 0.
.TP 8
A       (input) COMPLEX*16 array, dimension
(LDA,M) if SIDE = \(aqL\(aq,
(LDA,N) if SIDE = \(aqR\(aq
The i-th row must contain the vector which defines the
elementary reflector H(i), for i = 1,2,...,k, as returned by
ZTZRZF in the last k rows of its array argument A.
A is modified by the routine but restored on exit.
.TP 8
LDA     (input) INTEGER
The leading dimension of the array A. LDA >= max(1,K).
.TP 8
TAU     (input) COMPLEX*16 array, dimension (K)
TAU(i) must contain the scalar factor of the elementary
reflector H(i), as returned by ZTZRZF.
.TP 8
C       (input/output) COMPLEX*16 array, dimension (LDC,N)
On entry, the M-by-N matrix C.
On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
.TP 8
LDC     (input) INTEGER
The leading dimension of the array C. LDC >= max(1,M).
.TP 8
WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
.TP 8
LWORK   (input) INTEGER
The dimension of the array WORK.
If SIDE = \(aqL\(aq, LWORK >= max(1,N);
if SIDE = \(aqR\(aq, LWORK >= max(1,M).
For optimum performance LWORK >= N*NB if SIDE = \(aqL\(aq, and
LWORK >= M*NB if SIDE = \(aqR\(aq, where NB is the optimal
blocksize.

If LWORK = -1, then a workspace query is assumed; the routine
only calculates the optimal size of the WORK array, returns
this value as the first entry of the WORK array, and no error
message related to LWORK is issued by XERBLA.
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
.SH FURTHER DETAILS
Based on contributions by
.br
  A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA

