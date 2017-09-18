
/* TODO
 * 1) try to store data for complex2real fft, and perform fft on host (less data to transfer)
 * 2) use static polar optimization (from constant memory ?)
 * 3) use a for loop in m-direction to re-use threads at larger m's.
 * 4) use double2 (double4 ?) vector types to avoid the use of shared memory to write mangled data.
 * 5) multi-stream / multi-gpu computation ?
 * 6) allow several variants, which may change occupancy for large sizes.
 */
 
// NOTE variables gridDim.x, blockIdx.x, blockDim.x, threadIdx.x, and warpSize are defined in device functions

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include "sht_private.h"

#include <cufft.h>

// 256 for scalar SH_to_spat seems best on kepler.
#define THREADS_PER_BLOCK 256
// number of latitudes per thread:
#define NWAY 1
// number of l blocks to pre-compute
#define LSPAN 16
// the warp size is always 32 on cuda devices (up to Pascal at least)
#define WARPSZE 32
// maximum streams:
#define MAX_STRM 2

// adjustment for cuda
#undef SHT_L_RESCALE_FLY
#undef SHT_ACCURACY
#define SHT_L_RESCALE_FLY 1800
#define SHT_ACCURACY 1.0e-40

cudaStream_t strm[MAX_STRM];

extern "C"
void* shtns_malloc(size_t size) {
    void* ptr = NULL;
    cudaMallocHost(&ptr, size);		// allocate pinned memory (for faster transfers !)
    return ptr;
}

extern "C"
void shtns_free(void* p) {
    cudaFreeHost(p);
}

/// On KEPLER, This kernel is fastest with THREADS_PER_BLOCK=256 and NWAY=1
template<int S> __global__ void
leg_m0(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2)
{
    // im = 0
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int j = threadIdx.x;

    __shared__ double ak[THREADS_PER_BLOCK];
    __shared__ double qk[THREADS_PER_BLOCK/2];
    ak[j] = al[j];
    if ((j <= llim)&&(j<THREADS_PER_BLOCK/2)) qk[j] = ql[2*j];
    __syncthreads();

    int l = 0;
    int k = 0;	int kq = 0;
    double cost[NWAY];
    double y0[NWAY];    double y1[NWAY];
    double re[NWAY];    double ro[NWAY];

    for (int i=0; i<NWAY; i++) {
	cost[i] = (it+i<nlat_2) ? ct[it+i] : 0.0;
	y0[i] = ak[0];
	if (S==1) y0[i] *= rsqrt(1.0 - cost[i]*cost[i]);	// for vectors, divide by sin(theta)
    }
    for (int i=0; i<NWAY; i++) {
	re[i] = y0[i] * qk[0];
	y1[i] = y0[i] * ak[1] * cost[i];
    }
    for (int i=0; i<NWAY; i++) {
	ro[i] = y1[i] * qk[1];
    }
    al+=2;    l+=2;	k+=2;	kq+=2;
    while(l<llim) {
	if (k+6 >= THREADS_PER_BLOCK) {
	    __syncthreads();
	    ak[j] = al[j];
	    if ((j <= llim)&&(j<THREADS_PER_BLOCK/2)) qk[j] = ql[2*(l+j)];
	    k=0;	kq=0;
    	    __syncthreads();
	}
	for (int i=0; i<NWAY; i++)	y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
	for (int i=0; i<NWAY; i++)	re[i] += y0[i] * qk[kq];
	for (int i=0; i<NWAY; i++)	y1[i]  = ak[k+3]*cost[i]*y0[i] + ak[k+2]*y1[i];
	for (int i=0; i<NWAY; i++)	ro[i] += y1[i] * qk[kq+1];
	al+=4;	l+=2;	k+=4;	kq+=2;
    }
    if (l==llim) {
	for (int i=0; i<NWAY; i++)	y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
	for (int i=0; i<NWAY; i++)	re[i] += y0[i] * qk[kq];
    }

    for (int i=0; i<NWAY; i++) {
	if (it+i < nlat_2) {
	    q[it+i] = re[i]+ro[i];
	    q[nlat_2*2-1-(it+i)] = re[i]-ro[i];
	}
    }
/*
    if (it < nlat_2) {
        int l = 0;
	double cost = ct[it];
        double y0 = al[0];
        double re = y0 * ql[0];
        double y1 = y0 * al[1] * cost;
        double ro = y1 * ql[1];
        al+=2;    l+=2;
        while(l<llim) {
            y0  = al[1]*(cost*y1) + al[0]*y0;
            re += y0 * ql[l];
            y1  = al[3]*(cost*y0) + al[2]*y1;
            ro += y1 * ql[l+1];
            al+=4;	l+=2;
        }
        if (l==llim) {
            y0  = al[1]*cost*y1 + al[0]*y0;
            re += y0 * ql[l];
        }

        q[it] = re+ro;
        q[nlat_2*2-1-it] = re-ro;
    }
    */
}

__inline__ __device__
void warp_reduce_add_4(double& re, double& ro, double& ie, double& io) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
    re += __shfl_down(re, offset);
    ro += __shfl_down(ro, offset);
    ie += __shfl_down(ie, offset);
    io += __shfl_down(io, offset);
  }
}

__inline__ __device__
void warp_reduce_add_2(double& ev, double& od) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
    ev += __shfl_down(ev, offset);
    od += __shfl_down(od, offset);
  }
}

__inline__ __device__
void warp_reduce_add(double& ev) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
    ev += __shfl_down(ev, offset);
  }
}

#if (__CUDACC_VER_MAJOR__ < 8) || ( defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600 )
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
	old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif


/// THREADS_PER_BLOCK/LSPAN must be a power of 2 and <= WARPSZE
/// LSPAN must be a multiple of 2.
template<int S> __global__ void
ileg_m0(const double *al, const double *ct, const double *q, double *ql, const int llim, const int nlat_2)
{
    // im = 0
    const int it = (blockDim.x * blockIdx.x + threadIdx.x)*NWAY;
    const int j = threadIdx.x;

    __shared__ double ak[2*LSPAN+2];	// cache
    __shared__ double yl[LSPAN][THREADS_PER_BLOCK+1];	// padding to avoid bank conflicts
    __shared__ double reo[2][THREADS_PER_BLOCK+1];	// padding to avoid bank conflicts
    double cost, y0, y1;

    y0 = (it < nlat_2) ? q[it] : 0.0;		// north
    y1 = (it < nlat_2) ? q[nlat_2*2-1 - it] : 0.0;	// south
    reo[0][j] = y0+y1;				// even
    reo[1][j] = y0-y1;		// odd

    if (j < 2*LSPAN+2) ak[j] = al[j];
    #if THREADS_PER_BLOCK > WARPSZE
    __syncthreads();
    #endif

    int l = 0;
    if (it < nlat_2) {
	y0 = ct[it + nlat_2];		// weights are stored just after ct.
	cost = ct[it];
    } else {
	y0 = 0.0;	cost = 0.0;
    }
    if (S==1) y0 *= rsqrt(1.0 - cost*cost);	// for vectors, divide by sin(theta)
    y0 *= ak[0];
    y1 = y0 * ak[1] * cost;
    yl[0][j] = y0;
    yl[1][j] = y1;
    al+=2;
    while (l < llim) {
	for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
	    yl[k][j]     = y0;
	    y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
	    yl[k+1][j] = y1;
	    y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
	    al += 4;
	}

	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	double qll = 0.0;	// accumulator
	// now re-assign each thread an l (transpose)
	const int ll = j / (THREADS_PER_BLOCK/LSPAN);
	for (int i=0; i<THREADS_PER_BLOCK; i+= THREADS_PER_BLOCK/LSPAN) {
	    int it = j % (THREADS_PER_BLOCK/LSPAN) + i;
	    qll += reo[ll&1][it] * yl[ll][it];
	}
	// reduce_add within same l must be in same warp too:
	#if THREADS_PER_BLOCK/LSPAN > WARPSZE
	    #error "THREADS_PER_BLOCK/LSPAN > WARPSZE"
	#endif
	for (int ofs = THREADS_PER_BLOCK/(LSPAN*2); ofs > 0; ofs>>=1) {
	    qll += __shfl_down(qll, ofs, THREADS_PER_BLOCK/LSPAN);
	}
	if ( ((j % (THREADS_PER_BLOCK/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
	    if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
		ql[2*(l+ll)] = qll;
	    } else {
		atomicAdd(ql+2*(l+ll), qll);		// VERY slow atomic add on Kepler.
	    }
	}
	if (j<2*LSPAN) ak[j+2] = al[j];
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	l+=LSPAN;
    }
}


/** \internal convert from vector SH to scalar SH
    Vlm =  st*d(Slm)/dtheta + I*m*Tlm
    Wlm = -st*d(Tlm)/dtheta + I*m*Slm
*/
__global__ void
sphtor2scal_gpu(const double *mx, const double *slm, const double *tlm, double *vlm, double *wlm, const int llim, const int lmax, const int mres)
{
    const int ll = (blockDim.x-4) * blockIdx.x + threadIdx.x - 2;		// = 2*l + ((imag) ? 1 : 0)
    const int j = threadIdx.x;
    const int im = blockIdx.y;

    __shared__ double sl[THREADS_PER_BLOCK];
    __shared__ double tl[THREADS_PER_BLOCK];
    __shared__ double M[THREADS_PER_BLOCK];

    const int m = __mul24(im, mres);
    //int ofs = im*(2*(lmax+1) -m + mres);
//    const int xchg = 1 - 2*(ll&1);	// +1 for real and -1 for imag
    const int xchg = ll - (ll^1);	// -1 for real and +1 for imag
    const int ofs = __mul24( im, (((lmax+1)<<1) -m + mres) ) + ll;

    if ( (ll >= 0) && (ll < 2*(llim+1-m)) ) {
	M[j] = mx[ofs];
	sl[j] = slm[ofs];
	tl[j] = tlm[ofs];
    } else {
	M[j] = 0.0;
	sl[j] = 0.0;
	tl[j] = 0.0;
    }
    const double mimag = im * mres * (ll - (ll^1));

    __syncthreads();

//    if ((j>=2) && (j<THREADS_PER_BLOCK-2) && (ll < 2*(llim+2-m))) {
    if ((j<THREADS_PER_BLOCK-4) && (ll < 2*(llim+1-m))) {
	double ml = M[2*(j>>1)+1];
	double mu = M[2*(j>>1)+2];
	double v = mimag*tl[(j+2)^1]  +  (ml*sl[j] + mu*sl[j+4]);
	double w = mimag*sl[(j+2)^1]  -  (ml*tl[j] + mu*tl[j+4]);
	vlm[ofs+2*im+2] = v;
	wlm[ofs+2*im+2] = w;
    }
}

/** \internal convert from 2 scalar SH to vector SH
    Slm = - (I*m*Wlm + MX*Vlm) / (l*(l+1))
    Tlm = - (I*m*Vlm - MX*Wlm) / (l*(l+1))
**/
__global__ void
scal2sphtor_gpu(const double *mx, const double *vlm, const double *wlm, double *slm, double *tlm, const int llim, const int lmax, const int mres)
{
    const int ll = (blockDim.x-4) * blockIdx.x + threadIdx.x - 2;		// = 2*l + ((imag) ? 1 : 0)
    const int j = threadIdx.x;
    const int im = blockIdx.y;

    __shared__ double vl[THREADS_PER_BLOCK];
    __shared__ double wl[THREADS_PER_BLOCK];
    __shared__ double M[THREADS_PER_BLOCK];

    const int m = __mul24(im, mres);
    //const int xchg = 1 - 2*(j&1);	// +1 for real and -1 for imag
    const int xchg = (j^1) - j;		// +1 for real and -1 for imag
    int ofs = im*(2*(lmax+1) -m + mres)  + ll;
    //int ofs = __mul24( im, (((lmax+1)<<1) -m + mres) )  + ll;

    if ( (ll >= 0) && (ll < 2*(llim+1-m)) ) {
	M[j] = mx[ofs];
    } else M[j] = 0.0;

    if ( (ll >= 0) && (ll < 2*(llim+2-m)) ) {
	vl[j] = vlm[ofs+2*im];
	wl[j] = wlm[ofs+2*im];
    } else {
	vl[j] = 0.0;
	wl[j] = 0.0;
    }

    int ell = (ll>>1) + m + 1;		// +1 because we shift below

    __syncthreads();

//    if ((j>=2) && (j<THREADS_PER_BLOCK-2) && (ll < 2*(llim+1-m))) {
    if (j<THREADS_PER_BLOCK-4) {
	if ((ell <= llim) && (ell>0)) {
	    const double mimag = im * xchg*mres;
	    double ll_1 = 1.0 / __mul24(ell,ell+1);
	    double ml = M[2*(j>>1)+1];
	    double mu = M[2*(j>>1)+2];
	    double s = mimag*wl[(j+2)^1]  -  (ml*vl[j] + mu*vl[j+4]);
	    double t = mimag*vl[(j+2)^1]  +  (ml*wl[j] + mu*wl[j+4]);
	    slm[ofs+2] = s * ll_1;
	    tlm[ofs+2] = t * ll_1;
	} else if (ell <= lmax) {	// fill with zeros up to lmax (and l=0 too).
	    slm[ofs+2] = 0.0;
	    tlm[ofs+2] = 0.0;
	}
    }
}


/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim MUST BE <= 1800
/// S can only be 0 (for scalar) or 1 (for spin 1 / vector)
template<int S> __global__ void
leg_m_lowllim(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi)
{
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int im = blockIdx.y;
    const int j = threadIdx.x;
    const int m_inc = 2*nlat_2;
    const int k_inc = 1;

    // two arrays in shared memory of size blockDim.x :
    extern __shared__ double ak[];
    double* const qk = ak + blockDim.x;

    const double cost = (it < nlat_2) ? ct[it] : 0.0;

    if (im==0) {
	ak[j] = al[j];
	if (j<2*(llim+1)) qk[j] = ql[j];
	__syncthreads();
	int l = 0;
	int ka = 0;	int kq = 0;
	double y0 = ak[0];
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);	// for vectors, divide by sin(theta)
	double re = y0 * qk[0];
	double y1 = y0 * ak[1] * cost;
	double ro = y1 * qk[2];
	al+=2;    l+=2;		ka+=2;	kq+=2;
	while(l<llim) {
	    if (ka+6 >= blockDim.x) {
		__syncthreads();  
		ak[j] = al[j];
		qk[j] = ql[2*l+j];
		ka=0;	kq=0;
		__syncthreads();
	    }
	    y0  = ak[ka+1]*cost*y1 + ak[ka]*y0;
	    re += y0 * qk[2*kq];
	    y1  = ak[ka+3]*cost*y0 + ak[ka+2]*y1;
	    ro += y1 * qk[2*kq+2];
	    al+=4;	l+=2;	  ka+=4;    kq+=2;
	}
	if (l==llim) {
	    y0  = ak[ka+1]*cost*y1 + ak[ka]*y0;
	    re += y0 * qk[2*kq];
	}
	if (it<nlat_2) {
	    // store mangled for complex fft
	    q[it*k_inc] = re+ro;
	    q[(nlat_2*2-1-it)*k_inc] = re-ro;
	}
    } else { 	// m>0
	double rer,ror, rei, roi, y0, y1;
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

	y1 = sqrt(1.0 - cost*cost);		// y1 = sin(theta)
	ak[j] = al[j];
	qk[j] = ql[2*m+j];

	ror = 0.0;		roi = 0.0;
	rer = 0.0;		rei = 0.0;
	y0 = 1.0;
	l = m - S;
	do {		// sin(theta)^(m-S)
	    if (l&1) y0 *= y1;
	    y1 *= y1;
	} while(l >>= 1);
	
	__syncthreads();
	y0 *= ak[0];
	y1 = ak[1]*y0*cost;

	int ka = 2;
	l=m;		al+=2;
	int kq = 0;

	while (l<llim) {	// compute even and odd parts
	    if (2*kq+6 > blockDim.x) {
		__syncthreads();
		ak[j] = al[j];
		qk[j] = ql[2*l+j];
		ka=0;	kq=0;
		__syncthreads();
	    }
	    rer += y0 * qk[2*kq];	// real
	    rei += y0 * qk[2*kq+1];	// imag
	    y0 = ak[ka+1]*(cost*y1) + ak[ka]*y0;
	    ror += y1 * qk[2*kq+2];	// real
	    roi += y1 * qk[2*kq+3];	// imag
	    y1 = ak[ka+3]*(cost*y0) + ak[ka+2]*y1;
	    l+=2;	al+=4;	 ka+=4;	  kq+=2;
	}
	if (l==llim) {
	    rer += y0 * qk[2*kq];
	    rei += y0 * qk[2*kq+1];
	}

	/// store mangled for complex fft
	double nr = rer+ror;
	double sr = rer-ror;
	const double sgn = 1 - 2*(j&1);
	rei = __shfl_xor(rei, 1);
	roi = __shfl_xor(roi, 1);
	double nix = sgn*(rei+roi);
	double six = sgn*(rei-roi);
	if (it < nlat_2) {
	    q[im*m_inc + it*k_inc]                     = nr - nix;
	    q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
	    q[im*m_inc + (nlat_2*2-1-it)*k_inc]        = sr + six;
	    q[(nphi-im)*m_inc + (nlat_2*2-1-it)*k_inc] = sr - six;
	}
    }
}

/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim can be arbitrarily large (> 1800)
template<int S> __global__ void
leg_m_highllim(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi)
{
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int im = blockIdx.y;
    const int j = threadIdx.x;
    const int m_inc = 2*nlat_2;
    const int k_inc = 1;

    __shared__ double ak[THREADS_PER_BLOCK];	// cache
    __shared__ double qk[THREADS_PER_BLOCK];

    const double cost = (it < nlat_2) ? ct[it] : 0.0;

    if (im==0) {
	int l = 0;
	double y0 = al[0];
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);
	double re = y0 * ql[0];
	double y1 = y0 * al[1] * cost;
	double ro = y1 * ql[2];
	al+=2;    l+=2;
	while(l<llim) {
	    y0  = al[1]*(cost*y1) + al[0]*y0;
	    re += y0 * ql[2*l];
	    y1  = al[3]*(cost*y0) + al[2]*y1;
	    ro += y1 * ql[2*l+2];
	    al+=4;	l+=2;
	}
	if (l==llim) {
	    y0  = al[1]*cost*y1 + al[0]*y0;
	    re += y0 * ql[2*l];
	}
	if (it < nlat_2) {
	    // store mangled for complex fft
	    q[it*k_inc] = re+ro;
	    q[(nlat_2*2-1-it)*k_inc] = re-ro;
	}
    } else { 	// m>0
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*(l + S*im);
	double rer,ror, rei,roi, y0, y1;
	ror = 0.0;	roi = 0.0;
	rer = 0.0;	rei = 0.0;
	y1 = sqrt(1.0 - cost*cost);	// sin(theta)
	if (__any(m - llim*y1 <= max(50,llim/200))) {		// polar optimization (see Reinecke 2013), avoiding warp divergence
	    y0 = 1.0;	// y0
	    l = m - S;
	    int ny = 0;
	    int nsint = 0;
	    do {		// sin(theta)^(m-S)		(use rescaling to avoid underflow)
		if (l&1) {
		    y0 *= y1;
		    ny += nsint;
		    if (__any(y0 < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR))) {		// avoid warp divergence
			ny--;
			y0 *= SHT_SCALE_FACTOR;
		    }
		}
		y1 *= y1;
		nsint += nsint;
		if (__any(y1 < 1.0/SHT_SCALE_FACTOR)) {		// avoid warp divergence
		    nsint--;
		    y1 *= SHT_SCALE_FACTOR;
		}
	    } while(l >>= 1);
	    y0 *= al[0];
	    y1 = 0.0;
//	    y1 = al[1]*y0*cost;

	    l=m;	int ka = WARPSZE;
	    const int ofs = j & 0xFFE0;

	    while ( __all(ny<0) && (l<llim) ) {
		if (ka+4 >= WARPSZE) {
		    ak[j] = al[(j&31)];
		    ka=0;
		}
		y1 = ak[ka+1+ofs]*cost*y0 + ak[ka+ofs]*y1;
		y0 = ak[ka+3+ofs]*cost*y1 + ak[ka+2+ofs]*y0;
		l+=2;	al+=4;	ka+=4;
		if (fabs(y1) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0)
		{	// rescale when value is significant
		    ++ny;
		    y0 *= 1.0/SHT_SCALE_FACTOR;
		    y1 *= 1.0/SHT_SCALE_FACTOR;
		}
	    }
	    
	    ka = WARPSZE;
	    while (l<llim) {
		if (ka+4 >= WARPSZE) {		// cache coefficients
		    ak[j] = al[(j&31)];
		    qk[j] = ql[2*l+(j&31)];
		    ka = 0;
		}
		y1 = ak[ka+1+ofs]*cost*y0 + ak[ka+ofs]*y1;
		if (ny==0) {
		    rer += y0 * qk[ka+ofs];	// real
		    rei += y0 * qk[ka+1+ofs];	// imag
		    ror += y1 * qk[ka+2+ofs];	// real
		    roi += y1 * qk[ka+3+ofs];	// imag
		}
		else if (fabs(y0) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0)
		{	// rescale when value is significant
		    ++ny;
		    y0 *= 1.0/SHT_SCALE_FACTOR;
		    y1 *= 1.0/SHT_SCALE_FACTOR;
		}
		l+=2;	al+=4;
		y0 = ak[ka+3+ofs]*cost*y1 + ak[ka+2+ofs]*y0;
		ka+=4;
	    }
	    if ((l==llim) && (ny==0)) {
		rer += y0 * ql[2*l];
		rei += y0 * ql[2*l+1];
	    }
	}

	/// store mangled for complex fft
	double nr = rer+ror;
	double sr = rer-ror;
	const double sgn = 1 - 2*(j&1);
	rei = __shfl_xor(rei, 1);
	roi = __shfl_xor(roi, 1);
	double nix = sgn*(rei+roi);
	double six = sgn*(rei-roi);
	if (it < nlat_2) {
	    q[im*m_inc + it*k_inc]                     = nr - nix;
	    q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
	    q[im*m_inc + (nlat_2*2-1-it)*k_inc]        = sr + six;
	    q[(nphi-im)*m_inc + (nlat_2*2-1-it)*k_inc] = sr - six;
	}
    }
}


template<int S> __global__ void
ileg_m_lowllim(const double *al, const double *ct, const double *q, double *ql, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi)
{
    const int it = (blockDim.x * blockIdx.x + threadIdx.x)*NWAY;
    const int j = threadIdx.x;
    const int im = blockIdx.y;
    const int m_inc = 2*nlat_2;
//    const int k_inc = 1;

    __shared__ double ak[2*LSPAN+2];	// cache
    __shared__ double yl[LSPAN*THREADS_PER_BLOCK];
    __shared__ double reo[4*THREADS_PER_BLOCK];
    const int l_inc = THREADS_PER_BLOCK;
    const double cost = (it < nlat_2) ? ct[it] : 0.0;
    double y0, y1;


    if (im == 0) {
	if (j < 2*LSPAN+2) ak[j] = al[j];
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	y0 = (it < nlat_2) ? q[it] : 0.0;		// north
	y1 = (it < nlat_2) ? q[nlat_2*2-1 - it] : 0.0;	// south
	reo[j] = y0+y1;				// even
	reo[THREADS_PER_BLOCK +j] = y0-y1;		// odd

	int l = 0;
	y0 = (it < nlat_2) ? ct[it + nlat_2] : 0.0;		// weights are stored just after ct.
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);
	y0 *= ak[0];
	y1 = y0 * ak[1] * cost;
	yl[j] = y0;
	yl[l_inc +j] = y1;
	al+=2;
	while (l <= llim) {
	    for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
		yl[k*l_inc +j]     = y0;
		y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
		yl[(k+1)*l_inc +j] = y1;
		y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
		al += 4;
	    }
	    #if THREADS_PER_BLOCK > WARPSZE
	    __syncthreads();
	    #endif
	    double qll = 0.0;	// accumulator
	    // now re-assign each thread an l (transpose)
	    const int ll = j / (THREADS_PER_BLOCK/LSPAN);
	    for (int i=0; i<THREADS_PER_BLOCK; i+= THREADS_PER_BLOCK/LSPAN) {
		int it = j % (THREADS_PER_BLOCK/LSPAN) + i;
		qll += reo[(ll&1)*THREADS_PER_BLOCK +it] * yl[ll*l_inc +it];
	    }
	    // reduce_add within same l must be in same warp too:
	    #if THREADS_PER_BLOCK/LSPAN > WARPSZE
		#error "THREADS_PER_BLOCK/LSPAN > WARPSZE"
	    #endif
	    for (int ofs = THREADS_PER_BLOCK/(LSPAN*2); ofs > 0; ofs>>=1) {
		qll += __shfl_down(qll, ofs, THREADS_PER_BLOCK/LSPAN);
	    }
	    if ( ((j % (THREADS_PER_BLOCK/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
		if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
		    ql[2*(l+ll)] = qll;
		} else {
		    atomicAdd(ql+2*(l+ll), qll);		// VERY slow atomic add on Kepler.
		}
	    }
	    if (j<2*LSPAN) ak[j+2] = al[j];
	    #if THREADS_PER_BLOCK > WARPSZE
	    __syncthreads();
	    #endif
	    l+=LSPAN;
	}
    } else {	// im > 0
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

	if (j < 2*LSPAN+2) ak[j] = al[j];
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	const double sgn = 2*(j&1) - 1;	// -/+
	y0    = (it < nlat_2) ? q[im*m_inc + it] : 0.0;		// north imag (ani)
	double qer    = (it < nlat_2) ? q[(nphi-im)*m_inc + it] : 0.0;	// north real (an)
	y1    = (it < nlat_2) ? q[im*m_inc + nlat_2*2-1-it] : 0.0;	// south imag (asi)
	double qor    = (it < nlat_2) ? q[(nphi-im)*m_inc + nlat_2*2-1-it] : 0.0;	// south real (as)
	double qei = y0-qer;		qer += y0;		// ani = -qei[lane+1],   bni = qei[lane-1]
	double qoi = y1-qor;		qor += y1;		// bsi = -qoi[lane-1],   asi = qoi[lane+1];
	y0 = __shfl_xor(qei, 1);	// exchange between adjacent lanes.
	y1 = __shfl_xor(qoi, 1);
	reo[j] 			    = qer + qor;	// rer
	reo[THREADS_PER_BLOCK +j]   = qer - qor;	// ror
	reo[2*THREADS_PER_BLOCK +j] = sgn*(y0 - y1);	// rei
	reo[3*THREADS_PER_BLOCK +j] = sgn*(y0 + y1);	// roi
    	
	y1 = sqrt(1.0 - cost*cost);	// sin(theta)

	    y0 = 0.5 * ak[0];	// y0
	    l = m - S;
	    do {		// sin(theta)^(m-S)
		if (l&1) y0 *= y1;
		y1 *= y1;
	    } while(l >>= 1);
	    if (it < nlat_2)     y0 *= ct[it + nlat_2];		// include quadrature weights.
	    y1 = ak[1]*y0*cost;

	    l=m;		al+=2;
	    while (l <= llim) {
		for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
		    yl[k*l_inc +j]     = y0;
		    y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
		    yl[(k+1)*l_inc +j] = y1;
		    y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
		    al += 4;
		}

		#if THREADS_PER_BLOCK > WARPSZE
		__syncthreads();
		#endif
		double qlri = 0.0;	// accumulator
		// now re-assign each thread an l (transpose)
		const int ll = j / (THREADS_PER_BLOCK/LSPAN);
		const int ri = j / (THREADS_PER_BLOCK/(2*LSPAN)) % 2;	// real (0) or imag (1)
		for (int i=0; i<THREADS_PER_BLOCK; i+= THREADS_PER_BLOCK/(2*LSPAN)) {
		    int it = j % (THREADS_PER_BLOCK/(2*LSPAN)) + i;
		    qlri += reo[((ll&1)+2*ri)*THREADS_PER_BLOCK +it]   * yl[ll*l_inc +it];
		}
		// reduce_add within same l must be in same warp too:
		#if THREADS_PER_BLOCK/(2*LSPAN) > WARPSZE
		    #error "THREADS_PER_BLOCK/(2*LSPAN) > WARPSZE"
		#endif
		for (int ofs = THREADS_PER_BLOCK/(LSPAN*4); ofs > 0; ofs>>=1) {
		    qlri += __shfl_down(qlri, ofs, THREADS_PER_BLOCK/(LSPAN*2));
		}
		if ( ((j % (THREADS_PER_BLOCK/(2*LSPAN))) == 0) && ((l+ll)<=llim) ) {	// write result
		    if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
			ql[2*(l+ll)+ri]   = qlri;
		    } else {
			atomicAdd(ql+2*(l+ll)+ri, qlri);		// VERY slow atomic add on Kepler.
		    }
		}
		if (j<2*LSPAN) ak[j+2] = al[j];
		#if THREADS_PER_BLOCK > WARPSZE
		__syncthreads();
		#endif
		l+=LSPAN;
	    }
    }
}


template<int S> __global__ void
ileg_m_highllim(const double *al, const double *ct, const double *q, double *ql, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi)
{
    const int it = (blockDim.x * blockIdx.x + threadIdx.x)*NWAY;
    const int j = threadIdx.x;
    const int im = blockIdx.y;
    const int m_inc = 2*nlat_2;
//    const int k_inc = 1;

    __shared__ double ak[2*LSPAN+2];	// cache
    __shared__ double yl[LSPAN*THREADS_PER_BLOCK];
    __shared__ double reo[4*THREADS_PER_BLOCK];
    const int l_inc = THREADS_PER_BLOCK;
    const double cost = (it < nlat_2) ? ct[it] : 0.0;
    double y0, y1;


    if (im == 0) {
	if (j < 2*LSPAN+2) ak[j] = al[j];
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	y0 = (it < nlat_2) ? q[it] : 0.0;		// north
	y1 = (it < nlat_2) ? q[nlat_2*2-1 - it] : 0.0;	// south
	reo[j] = y0+y1;				// even
	reo[THREADS_PER_BLOCK +j] = y0-y1;		// odd

	int l = 0;
	y0 = (it < nlat_2) ? ct[it + nlat_2] : 0.0;		// weights are stored just after ct.
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);
	y0 *= ak[0];
	y1 = y0 * ak[1] * cost;
	yl[j] = y0;
	yl[l_inc +j] = y1;
	al+=2;
	while (l <= llim) {
	    for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
		yl[k*l_inc +j]     = y0;
		y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
		yl[(k+1)*l_inc +j] = y1;
		y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
		al += 4;
	    }
	    #if THREADS_PER_BLOCK > WARPSZE
	    __syncthreads();
	    #endif
	    double qll = 0.0;	// accumulator
	    // now re-assign each thread an l (transpose)
	    const int ll = j / (THREADS_PER_BLOCK/LSPAN);
	    for (int i=0; i<THREADS_PER_BLOCK; i+= THREADS_PER_BLOCK/LSPAN) {
		int it = j % (THREADS_PER_BLOCK/LSPAN) + i;
		qll += reo[(ll&1)*THREADS_PER_BLOCK +it] * yl[ll*l_inc +it];
	    }
	    // reduce_add within same l must be in same warp too:
	    #if THREADS_PER_BLOCK/LSPAN > WARPSZE
		#error "THREADS_PER_BLOCK/LSPAN > WARPSZE"
	    #endif
	    for (int ofs = THREADS_PER_BLOCK/(LSPAN*2); ofs > 0; ofs>>=1) {
		qll += __shfl_down(qll, ofs, THREADS_PER_BLOCK/LSPAN);
	    }
	    if ( ((j % (THREADS_PER_BLOCK/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
		if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
		    ql[2*(l+ll)] = qll;
		} else {
		    atomicAdd(ql+2*(l+ll), qll);		// VERY slow atomic add on Kepler.
		}
	    }
	    if (j<2*LSPAN) ak[j+2] = al[j];
	    #if THREADS_PER_BLOCK > WARPSZE
	    __syncthreads();
	    #endif
	    l+=LSPAN;
	}
    } else {	// im > 0
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

	if (j < 2*LSPAN+2) ak[j] = al[j];
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	const double sgn = 2*(j&1) - 1;	// -/+
	y0    = (it < nlat_2) ? q[im*m_inc + it] : 0.0;		// north imag (ani)
	double qer    = (it < nlat_2) ? q[(nphi-im)*m_inc + it] : 0.0;	// north real (an)
	y1    = (it < nlat_2) ? q[im*m_inc + nlat_2*2-1-it] : 0.0;	// south imag (asi)
	double qor    = (it < nlat_2) ? q[(nphi-im)*m_inc + nlat_2*2-1-it] : 0.0;	// south real (as)
	double qei = y0-qer;		qer += y0;		// ani = -qei[lane+1],   bni = qei[lane-1]
	double qoi = y1-qor;		qor += y1;		// bsi = -qoi[lane-1],   asi = qoi[lane+1];
	y0 = __shfl_xor(qei, 1);	// exchange between adjacent lanes.
	y1 = __shfl_xor(qoi, 1);
	reo[j] 			    = qer + qor;	// rer
	reo[THREADS_PER_BLOCK +j]   = qer - qor;	// ror
	reo[2*THREADS_PER_BLOCK +j] = sgn*(y0 - y1);	// rei
	reo[3*THREADS_PER_BLOCK +j] = sgn*(y0 + y1);	// roi
    	
	y1 = sqrt(1.0 - cost*cost);	// sin(theta)

	    y0 = 0.5;	// y0
	    l = m - S;
	    int ny = 0;
	    int nsint = 0;
	    do {		// sin(theta)^(m-S)		(use rescaling to avoid underflow)
		if (l&1) {
		    y0 *= y1;
		    ny += nsint;
		    // the use of __any leads to wrong results. On KEPLER it is also slower.
//		    if (__any(y0 < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR))) {		// avoid warp divergence
		    if (y0 < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
			ny--;
			y0 *= SHT_SCALE_FACTOR;
		    }
		}
		y1 *= y1;
		nsint += nsint;
//		if (__any(y1 < 1.0/SHT_SCALE_FACTOR)) {		// avoid warp divergence
		if (y1 < 1.0/SHT_SCALE_FACTOR) {
		    nsint--;
		    y1 *= SHT_SCALE_FACTOR;
		}
	    } while(l >>= 1);
	    y0 *= ak[0];
	    if (it < nlat_2)     y0 *= ct[it + nlat_2];		// include quadrature weights.
	    y1 = ak[1]*y0*cost;


	    l=m;		al+=2;
	    while (l <= llim) {
		for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
		    yl[k*l_inc +j]     = (ny==0) ? y0 : 0.0;
		    y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
		    yl[(k+1)*l_inc +j] = (ny==0) ? y1 : 0.0;
		    y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
		    if (ny<0) {
//			if (__any(fabs(y0) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0))
			if (fabs(y0) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0)
			{	// rescale when value is significant
			    ++ny;
			    y0 *= 1.0/SHT_SCALE_FACTOR;
			    y1 *= 1.0/SHT_SCALE_FACTOR;
			}
		    }
		    al += 4;
		}

		#if THREADS_PER_BLOCK > WARPSZE
		__syncthreads();
		#endif
		double qlri = 0.0;	// accumulator
		// now re-assign each thread an l (transpose)
		const int ll = j / (THREADS_PER_BLOCK/LSPAN);
		const int ri = j / (THREADS_PER_BLOCK/(2*LSPAN)) % 2;	// real (0) or imag (1)
		if (ll+l <= llim) {
		    for (int i=0; i<THREADS_PER_BLOCK; i+= THREADS_PER_BLOCK/(2*LSPAN)) {
			int it = j % (THREADS_PER_BLOCK/(2*LSPAN)) + i;
			qlri += reo[((ll&1)+2*ri)*THREADS_PER_BLOCK +it]   * yl[ll*l_inc +it];
		    }
		}
		// reduce_add within same l must be in same warp too:
		#if THREADS_PER_BLOCK/(2*LSPAN) > WARPSZE
		    #error "THREADS_PER_BLOCK/(2*LSPAN) > WARPSZE"
		#endif
		for (int ofs = THREADS_PER_BLOCK/(LSPAN*4); ofs > 0; ofs>>=1) {
		    qlri += __shfl_down(qlri, ofs, THREADS_PER_BLOCK/(LSPAN*2));
		}
		if ( ((j % (THREADS_PER_BLOCK/(2*LSPAN))) == 0) && ((l+ll)<=llim) ) {	// write result
		    if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
			ql[2*(l+ll)+ri]   = qlri;
		    } else {
			atomicAdd(ql+2*(l+ll)+ri, qlri);		// VERY slow atomic add on Kepler.
		    }
		}
		if (j<2*LSPAN) ak[j+2] = al[j];
		#if THREADS_PER_BLOCK > WARPSZE
		__syncthreads();
		#endif
		l+=LSPAN;
	    }
    }
}





extern "C"
int cushtns_init_gpu(shtns_cfg shtns)
{
    cudaError_t err = cudaSuccess;
    const long nlm = shtns->nlm;
    const long nlat_2 = shtns->nlat_2;

    double *d_alm = NULL;
    double *d_ct  = NULL;
    double *d_qlm = NULL;
    double *d_q   = NULL;
    double *d_mx_stdt = NULL;
    double *d_mx_van = NULL;
    int err_count = 0;
    int device_id = -1;

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, 0);
    if (err != cudaSuccess) return -1;
    #if SHT_VERBOSE > 0
    printf("  cuda GPU \"%s\" found (warp size = %d, compute capabilities = %d.%d).\n", prop.name, prop.warpSize, prop.major, prop.minor);
    #endif
    if (prop.warpSize != WARPSZE) return -1;		// failure, SHTns requires a warpSize of 32.
    if (prop.major < 3) return -1;			// failure, SHTns requires compute cap. >= 3 (warp shuffle instructions)

    // Allocate the device input vector alm
    err = cudaMalloc((void **)&d_alm, (2*nlm+THREADS_PER_BLOCK-1)*sizeof(double));	// allow some overflow.
    if (err != cudaSuccess) err_count ++;
    if (shtns->mx_stdt) {
	// Allocate the device matrix for d(sin(t))/dt
	err = cudaMalloc((void **)&d_mx_stdt, (2*nlm+THREADS_PER_BLOCK-1)*sizeof(double));
	if (err != cudaSuccess) err_count ++;
	// Same thing for analysis
	err = cudaMalloc((void **)&d_mx_van, (2*nlm+THREADS_PER_BLOCK-1)*sizeof(double));
	if (err != cudaSuccess) err_count ++;
    }
    // Allocate the device input vector cos(theta) and gauss weights
    err = cudaMalloc((void **)&d_ct, 2*nlat_2*sizeof(double));
    if (err != cudaSuccess) err_count ++;
    // Allocate the device work vector qlm
    err = cudaMalloc((void **)&d_qlm, 2*nlm*sizeof(double) * MAX_STRM);
    if (err != cudaSuccess) err_count ++;
    // Allocate the device work vector q
    err = cudaMalloc((void **)&d_q, shtns->nlat * shtns->nphi * sizeof(double) * MAX_STRM);
    if (err != cudaSuccess) err_count ++;

    if (err_count == 0) {
	err = cudaMemcpy(d_alm, shtns->alm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
	if (shtns->mx_stdt) {
	    err = cudaMemcpy(d_mx_stdt, shtns->mx_stdt, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess)  err_count ++;
	    err = cudaMemcpy(d_mx_van, shtns->mx_van, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	    if (err != cudaSuccess)  err_count ++;
	}
	err = cudaMemcpy(d_ct, shtns->ct, nlat_2*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
	err = cudaMemcpy(d_ct + nlat_2, shtns->wg, nlat_2*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
    }
    
    for (int i=0; i<MAX_STRM; i++) {
	cudaStreamCreate(&strm[i]);
    }

    /* cuFFT init */
    int nfft = shtns->nphi;
    if (nfft > 1) {
	// cufftPlanMany(cufftHandle *plan, int rank, int *n,   int *inembed, int istride, int idist,   int *onembed, int ostride, int odist,   cufftType type, int batch);
	cufftResult res;
	for (int i=0; i<MAX_STRM; i++) {
	    res = cufftPlanMany((cufftHandle*) &shtns->cufft_plan[i], 1, &nfft, &nfft, shtns->nlat_2, 1, &nfft, shtns->nlat_2, 1, CUFFT_Z2Z, shtns->nlat_2);
	    if (res != CUFFT_SUCCESS)  err_count ++;
	    cufftSetStream(shtns->cufft_plan[i], strm[i]);
	}
    }

    if (err_count != 0) {
	cudaFree(d_q);	cudaFree(d_qlm);  cudaFree(d_mx_van);	cudaFree(d_mx_stdt);  cudaFree(d_ct);  cudaFree(d_alm);
	return -1;	// fail
    }

    shtns->d_alm = d_alm;
    shtns->d_ct  = d_ct;
    shtns->d_q   = d_q;
    shtns->d_qlm = d_qlm;
    shtns->d_mx_stdt = d_mx_stdt;
    shtns->d_mx_van = d_mx_van;
    cudaGetDevice(&device_id);
    return device_id;		// success, return device_id
}

extern "C"
void cushtns_release_gpu(shtns_cfg shtns)
{
    cufftDestroy(shtns->cufft_plan[0]);
    cufftDestroy(shtns->cufft_plan[1]);
    cudaStreamDestroy(strm[0]);
    cudaStreamDestroy(strm[1]);
    if (shtns->d_q) cudaFree(shtns->d_q);
    if (shtns->d_qlm) cudaFree(shtns->d_qlm);
    if (shtns->d_ct) cudaFree(shtns->d_ct);
    if (shtns->d_alm) cudaFree(shtns->d_alm);
    if (shtns->d_mx_stdt) cudaFree(shtns->d_mx_stdt);
    shtns->d_alm = 0;
}

/// \internal Enables parallel transforms on selected GPU device, if available. \see shtns_use_gpu 
extern "C"
int cushtns_use_gpu(int device_id)
{
    int count = 0;
    if (device_id >= 0) {
	cudaGetDeviceCount(&count);
	if (count > 0) {
	    device_id = device_id % count;
	    cudaSetDevice(device_id);
	    return device_id;
	}
    }
    return -1;		// disable gpu.
}


extern "C"
void SH_to_spat_gpu_hostfft(shtns_cfg shtns, cplx *Qlm, double *Vr, const long int llim)
{
    cudaError_t err = cudaSuccess;
    const int lmax = shtns->lmax;
    int mmax = shtns->mmax;
    const int mres = shtns->mres;
    const int nlm = shtns->nlm;
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat/2 + threadsPerBlock*NWAY - 1) / (threadsPerBlock*NWAY);
    double *d_qlm = shtns->d_qlm;
    double *d_q = shtns->d_q;
    if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
    if (mmax == 0) {
	double* Ql0;
	Ql0 = (double*) malloc((lmax+1)*sizeof(double));
	for (int l=0; l<=llim; l++) {
	    Ql0[l] = creal(Qlm[l]);
	}
	err = cudaMemcpy(d_qlm, Ql0, (llim+1)*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	leg_m0<0><<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2);
    } else {
	err = cudaMemcpy(d_qlm, Qlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<0><<<blocks, threads, 2*threadsPerBlock*sizeof(double)>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2, lmax,mres, nphi);
	} else {
	    leg_m_highllim<0><<<blocks, threads>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2, lmax,mres, nphi);
	}
	// padd missing m's with 0 (m>mmax)
	if (2*(mmax+1) <= nphi)
	    cudaMemset( d_q + (mmax+1)*nlat, 0, sizeof(double)*(nphi-2*mmax-1)*nlat );		// set to zero before fft
    }
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Failed to launch cuda kernel (error code %s)!\n", cudaGetErrorString(err));
    }

    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)  printf("failed copy back : %s\n", cudaGetErrorString(err));

    if (nphi > 1) {		// fft
  	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode == 0) {
		    fftw_execute_dft(shtns->ifftc, (cplx *) Vr, (cplx *) Vr);
		} else {		// split dft
		    printf("ERROR fft not supported\n");
		}
	}
    }
}

/// Perform SH transform on data that is already on the GPU. d_Qlm and d_Vr are pointers to GPU memory (obtained by cudaMalloc() for instance)
template<int S>
void cuda_SH_to_spat(shtns_cfg shtns, cplx* d_Qlm, double *d_Vr, const long int llim, int strm_idx = 0)
{
    const int lmax = shtns->lmax;
    int mmax = shtns->mmax;
    const int mres = shtns->mres;
    const int nlat_2 = shtns->nlat_2;
    const int nphi = shtns->nphi;
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;
    cudaStream_t stream = strm[strm_idx];

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat_2 + threadsPerBlock - 1) / threadsPerBlock;
    if (nphi == 1) {
	leg_m0<S><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2);
    } else {
	if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (mmax==0) {
	    leg_m0<S><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2);
	} else
	if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<S><<<blocks, threads, 2*threadsPerBlock*sizeof(double), stream>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mres, nphi);
	} else {
	    leg_m_highllim<S><<<blocks, threads, 0, stream>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mres, nphi);
	}
	// padd missing m's with 0 (m>mmax)
	if (2*(mmax+1) <= nphi)
	    cudaMemsetAsync( d_Vr + (mmax+1)*2*nlat_2, 0, sizeof(double)*(nphi-2*mmax-1)*2*nlat_2, stream );		// set to zero before fft
	cufftResult res;
	res = cufftExecZ2Z((cufftHandle) shtns->cufft_plan[strm_idx], (cufftDoubleComplex*) d_Vr, (cufftDoubleComplex*) d_Vr, CUFFT_INVERSE);
	if (res != CUFFT_SUCCESS) printf("cufft error %d\n", res);
    }
}

extern "C"
void cu_SH_to_spat(shtns_cfg shtns, cplx* d_Qlm, double *d_Vr, const long int llim)
{
    cuda_SH_to_spat<0>(shtns, d_Qlm, d_Vr, llim);
}

extern "C"
void cu_SHsphtor_to_spat(shtns_cfg shtns, cplx* d_Slm, cplx* d_Tlm, double* d_Vt, double* d_Vp, const long llim)
{
    double* d_vwlm;
    const int nlm = shtns->nlm + (shtns->mmax+1);	// we need one more mode per m.
    const long nlm_stride = ((2*nlm+WARPSZE-1)/WARPSZE) * WARPSZE;
    cudaError_t err = cudaSuccess;

    // we need temporary storage here ...
    err = cudaMalloc( (void **)&d_vwlm, (2*nlm_stride)*sizeof(double) );

    dim3 blocks((2*(shtns->lmax+2)+THREADS_PER_BLOCK-5)/(THREADS_PER_BLOCK-4), shtns->mmax+1);
    dim3 threads(THREADS_PER_BLOCK, 1);
    sphtor2scal_gpu <<<blocks, threads>>>
	(shtns->d_mx_stdt, (double*) d_Slm, (double*) d_Tlm, d_vwlm, d_vwlm+nlm_stride, llim, shtns->lmax, shtns->mres);
//    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("sphtor2scal_gpu error : %s!\n", cudaGetErrorString(err));	return; }

    // SHT on the GPU
    cuda_SH_to_spat<1>(shtns, (cplx*) d_vwlm, d_Vt, llim+1, 0);
    cuda_SH_to_spat<1>(shtns, (cplx*) (d_vwlm + nlm_stride), d_Vp, llim+1, 1);

    cudaFree(d_vwlm);
}

extern "C"
void cu_SHqst_to_spat(shtns_cfg shtns, cplx* d_Qlm, cplx* d_Slm, cplx* d_Tlm, double* d_Vr, double* d_Vt, double* d_Vp, const long llim)
{
    double* d_vwlm;
    const int nlm = shtns->nlm + (shtns->mmax+1);	// we need one more mode per m.
    const long nlm_stride = ((2*nlm+WARPSZE-1)/WARPSZE) * WARPSZE;
    cudaError_t err = cudaSuccess;

    cuda_SH_to_spat<0>(shtns, d_Qlm, d_Vr, llim, 0);	// scalar part on stream 0

    // we need temporary storage here ...
    err = cudaMalloc( (void **)&d_vwlm, (2*nlm_stride)*sizeof(double) );

    dim3 blocks((2*(shtns->lmax+2)+THREADS_PER_BLOCK-5)/(THREADS_PER_BLOCK-4), shtns->mmax+1);
    dim3 threads(THREADS_PER_BLOCK, 1);
    sphtor2scal_gpu <<<blocks, threads, 0, strm[1]>>>
	(shtns->d_mx_stdt, (double*) d_Slm, (double*) d_Tlm, d_vwlm, d_vwlm+nlm_stride, llim, shtns->lmax, shtns->mres);
    cudaDeviceSynchronize();		// make stream 2 wait for stream 1
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("sphtor2scal_gpu error : %s!\n", cudaGetErrorString(err));	return; }

    // SHT on the GPU
    cuda_SH_to_spat<1>(shtns, (cplx*) d_vwlm, d_Vt, llim+1, 1);			// stream 1
    cuda_SH_to_spat<1>(shtns, (cplx*) (d_vwlm + nlm_stride), d_Vp, llim+1, 2);	// stream 2

    cudaFree(d_vwlm);
}


/// Perform SH transform on data that is already on the GPU. d_Qlm and d_Vr are pointers to GPU memory (obtained by cudaMalloc() for instance)
template<int S>
void cuda_spat_to_SH(shtns_cfg shtns, double *d_Vr, cplx* d_Qlm, const long int llim, int strm_idx = 0)
{
    const int lmax = shtns->lmax;
    int mmax = shtns->mmax;
    const int mres = shtns->mres;
    const int nlat_2 = shtns->nlat_2;
    const int nphi = shtns->nphi;
    const int nlm = shtns->nlm +S*(mmax+1);	// use more space for vector transform !!!
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;
    cudaStream_t stream = strm[strm_idx];

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat_2 + threadsPerBlock - 1) / (threadsPerBlock*NWAY);
    cudaMemsetAsync(d_Qlm, 0, sizeof(double)*2*nlm, stream);		// set to zero before we start.
    if (nphi == 1) {
	ileg_m0<S><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2);
    } else {
	cufftResult res;
	res = cufftExecZ2Z((cufftHandle) shtns->cufft_plan[strm_idx], (cufftDoubleComplex*) d_Vr, (cufftDoubleComplex*) d_Vr, CUFFT_INVERSE);
	if (res != CUFFT_SUCCESS) printf("cufft error %d\n", res);

	if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (mmax==0) {
	    ileg_m0<S><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2);
	} else
	if (llim <= SHT_L_RESCALE_FLY) {
	    ileg_m_lowllim<S><<<blocks, threads, 0, stream>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2, lmax,mres,nphi);
	} else {
	    ileg_m_highllim<S><<<blocks, threads, 0, stream>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2, lmax,mres,nphi);
	}
    }
}

extern "C"
void cu_spat_to_SH(shtns_cfg shtns, double *d_Vr, cplx* d_Qlm, const long int llim)
{
    cuda_spat_to_SH<0>(shtns, d_Vr, d_Qlm, llim);
}


extern "C"
void SH_to_spat_gpu(shtns_cfg shtns, cplx *Qlm, double *Vr, const long int llim)
{
    cudaError_t err = cudaSuccess;
    const int nlm = shtns->nlm;
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;

    double *d_qlm = shtns->d_qlm;
    double *d_q = shtns->d_q;

    // Allocate the device work vectors qlm and q
//    err = cudaMalloc((void **)&d_qlm, ((2*nlm +31 + nlat*nphi+31)/32)*32*sizeof(double));
//    d_q = d_qlm + ((2*nlm+31)/32)*32;

    // copy spectral data to GPU
    err = cudaMemcpy(d_qlm, Qlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("SH_to_spat_gpu failed copy qlm\n");	return; }

    // SHT on the GPU
    cuda_SH_to_spat<0>(shtns, (cplx*) d_qlm, d_q, llim);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("SH_to_spat_gpu CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spatial data
    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { printf("SH_to_spat_gpu failed copy back: %s\n", cudaGetErrorString(err));	return; }
    
//    cudaFree(d_qlm);
}


/** \internal convert from vector SH to scalar SH
    Vlm =  st*d(Slm)/dtheta + I*m*Tlm
    Wlm = -st*d(Tlm)/dtheta + I*m*Slm
**/
void sphtor2scal(shtns_cfg shtns, cplx* Slm, cplx* Tlm, cplx* Vlm, cplx* Wlm, const int llim)
{
    const int mmax = shtns->mmax;
    const int lmax = shtns->lmax;
    const int mres = shtns->mres;
    for (int im=0; im<=mmax; im++) {
	const int m = im*mres;
	long l = (im*(2*(lmax+1)-(m+mres)))>>1;
	double* mx = shtns->mx_stdt + 2*l;
	cplx* Sl = (cplx*) &Slm[l];	// virtual pointer for l=0 and im
	cplx* Tl = (cplx*) &Tlm[l];
	cplx* Vl = (cplx*) &Vlm[l+im];
	cplx* Wl = (cplx*) &Wlm[l+im];
	const double em = m;
	
	cplx sl = Sl[m];
	cplx tl = Tl[m];
	cplx vs = 0.0;
	cplx wt = 0.0;
	for (int l=m; l<=llim; l++) {
	    double mxu = mx[2*l];
	    double mxl = mx[2*l+1];	// mxl for next l
	    vs += I*em*tl;
	    wt += I*em*sl;
	    cplx vs1 = mxl*sl;		// vs for next l
	    cplx wt1 = -mxl*tl;		// wt for next l
	    if (l<llim) {
		sl = Sl[l+1];		// kept for next l
		tl = Tl[l+1];
		vs += mxu*sl;
		wt -= mxu*tl;
	    }
	    Vl[l] = vs;
	    Wl[l] = wt;
	    vs = vs1;		wt = wt1;
	}
	Vl[llim+1] = vs;
	Wl[llim+1] = wt;
    }
}


extern "C"
void SHsphtor_to_spat_gpu(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, const long int llim)
{
    cudaError_t err = cudaSuccess;
    const int nlm = shtns->nlm;
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;
    const int nlm2 = nlm + (shtns->mmax+1);	// one more data per m

    static double* d_vwlm = NULL;
    double* d_vtp;
//    static double* vw = NULL;

    const long nlm_stride = ((2*nlm2+WARPSZE-1)/WARPSZE) * WARPSZE;
    const long spat_stride = ((nlat*nphi+WARPSZE-1)/WARPSZE) * WARPSZE;

//    if (vw == NULL)
//	err = cudaMallocHost( (void**) &vw, ((nlm_stride > spat_stride) ? nlm_stride : spat_stride)*2*sizeof(double) );	// pinned buffer for transfer
    // Allocate the device work vectors
    if (d_vwlm == NULL)
    err = cudaMalloc( (void **)&d_vwlm, (4*nlm_stride + 2*spat_stride)*sizeof(double) );
    d_vtp = d_vwlm + 4*nlm_stride;

/*   // convert on cpu & transfer (via pinned mem) 
    sphtor2scal(shtns, Slm, Tlm, (cplx*) vw, (cplx*) (vw + nlm_stride), llim);		// convert & copy to pinned mem
    err = cudaMemcpy(d_vwlm, vw, 2*nlm_stride*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("failed copy vw\n");	return; }
*/
    // OR transfer and convert on gpu
    err = cudaMemcpy(d_vwlm + 2*nlm_stride, Slm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("memcpy 1 error : %s!\n", cudaGetErrorString(err));	return; }
    err = cudaMemcpy(d_vwlm + 3*nlm_stride, Tlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("memcpy 2 error : %s!\n", cudaGetErrorString(err));	return; }
    dim3 blocks((2*(shtns->lmax+2)+THREADS_PER_BLOCK-5)/(THREADS_PER_BLOCK-4), shtns->mmax+1);
    dim3 threads(THREADS_PER_BLOCK, 1);
    sphtor2scal_gpu <<<blocks, threads>>>
	(shtns->d_mx_stdt, d_vwlm+2*nlm_stride, d_vwlm+3*nlm_stride, d_vwlm, d_vwlm+nlm_stride, llim, shtns->lmax, shtns->mres);
//    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("sphtor2scal_gpu error : %s!\n", cudaGetErrorString(err));	return; }

    // SHT on the GPU
    cuda_SH_to_spat<1>(shtns, (cplx*) d_vwlm, d_vtp, llim+1, 0);
    cuda_SH_to_spat<1>(shtns, (cplx*) (d_vwlm + nlm_stride), d_vtp + spat_stride, llim+1, 1);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("SH_to_spat CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spatial data
    err = cudaMemcpy(Vt, d_vtp, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(Vp, d_vtp + spat_stride, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);

/*	// OR copy to pinned memory first
    err = cudaMemcpy(vw, d_vtp, 2*spat_stride*sizeof(double), cudaMemcpyDeviceToHost);
    memcpy(Vt, vw, nlat*nphi*sizeof(double));
    memcpy(Vp, vw + spat_stride, nlat*nphi*sizeof(double));
*/

//    cudaFree(d_vwlm);		d_vwlm = NULL;
//    cudaFreeHost(vw);
}



extern "C"
void spat_to_SH_gpu(shtns_cfg shtns, double *Vr, cplx *Qlm, const long int llim)
{
    cudaError_t err = cudaSuccess;
    const int nlm = shtns->nlm;
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;

    double *d_qlm = shtns->d_qlm;
    double *d_q = shtns->d_q;

    // copy spatial data to GPU
    err = cudaMemcpy(d_q, Vr, nlat*nphi*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("spat_to_SH_gpu failed copy q\n");	return; }

    // SHT on the GPU
    cu_spat_to_SH(shtns, d_q, (cplx*) d_qlm, llim);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("spat_to_SH_gpu CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spectral data
    err = cudaMemcpy(Qlm, d_qlm, 2*nlm*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { printf("spat_to_SH_gpu failed copy back\n");	return; }
}


extern "C"
void spat_to_SHsphtor_gpu(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, const long int llim)
{
    cudaError_t err = cudaSuccess;
    const int nlm = shtns->nlm;
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;
    const int nlm2 = nlm + (shtns->mmax+1);	// one more data per m

    static double* d_vwlm = NULL;
    double* d_vtp;
//    static double* vw = NULL;

    const long nlm_stride = ((2*nlm2+WARPSZE-1)/WARPSZE) * WARPSZE;
    const long spat_stride = ((nlat*nphi+WARPSZE-1)/WARPSZE) * WARPSZE;

//    if (vw == NULL)
//	err = cudaMallocHost( (void**) &vw, ((nlm_stride > spat_stride) ? nlm_stride : spat_stride)*2*sizeof(double) );	// pinned buffer for transfer
    // Allocate the device work vectors
    if (d_vwlm == NULL)
	 err = cudaMalloc( (void **)&d_vwlm, (4*nlm_stride + 2*spat_stride)*sizeof(double) );
    d_vtp = d_vwlm + 4*nlm_stride;

    // copy spatial data to gpu
    err = cudaMemcpy(d_vtp, Vt, nlat*nphi*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("memcpy 3 error : %s!\n", cudaGetErrorString(err));	return; }
    err = cudaMemcpy(d_vtp + spat_stride, Vp, nlat*nphi*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("memcpy 4 error : %s!\n", cudaGetErrorString(err));	return; }

    // SHT on the GPU
    cuda_spat_to_SH<1>(shtns, d_vtp, (cplx*) d_vwlm, llim+1, 0);
    cuda_spat_to_SH<1>(shtns, d_vtp + spat_stride, (cplx*) (d_vwlm + nlm_stride), llim+1, 0);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("spat_to_SHsphtor CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    cudaDeviceSynchronize();

    dim3 blocks((2*(shtns->lmax+2)+THREADS_PER_BLOCK-5)/(THREADS_PER_BLOCK-4), shtns->mmax+1);
    dim3 threads(THREADS_PER_BLOCK, 1);
    scal2sphtor_gpu <<<blocks, threads>>>
	(shtns->d_mx_van, d_vwlm, d_vwlm+nlm_stride, d_vwlm+2*nlm_stride, d_vwlm+3*nlm_stride, llim, shtns->lmax, shtns->mres);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("scal2sphtor_gpu error : %s!\n", cudaGetErrorString(err));	return; }

    cudaDeviceSynchronize();

    err = cudaMemcpy(Slm, d_vwlm + 2*nlm_stride, 2*nlm*sizeof(double), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(Tlm, d_vwlm + 3*nlm_stride, 2*nlm*sizeof(double), cudaMemcpyDeviceToHost);

//    cudaFree(d_vwlm);		d_vwlm = NULL;
//    cudaFreeHost(vw);
}



extern "C"
void SH_to_spat_many_gpu(shtns_cfg shtns, int howmany, cplx *Qlm, double *Vr, const long int llim)
{
    cudaError_t err = cudaSuccess;
    cudaEvent_t* event;
    const int nlm = shtns->nlm;
    const int lmax = shtns->lmax;
    int mmax = shtns->mmax;
    const int mres = shtns->mres;    
    const int nlat = shtns->nlat;
    const int nphi = shtns->nphi;
    const int nspat = shtns->nspat;
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;

    double* pinned;

    double *d_qlm = shtns->d_qlm;
    double *d_q = shtns->d_q;

    const int nstreams = (howmany < MAX_STRM) ? howmany : MAX_STRM;

    size_t dist = (2*nlm > nspat) ? 2*nlm : nspat;		// largest size between input or output.
    dist = ((dist + 31) >> 5) << 5;	// round to 32*8 = 256 bytes.
    cudaMallocHost(&pinned, sizeof(double)*howmany*dist);	// alloc pinned memory for fast transfers.

    event = (cudaEvent_t*) malloc(sizeof(cudaEvent_t) * howmany);

    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat/2 + threadsPerBlock*NWAY - 1) / (threadsPerBlock*NWAY);

    // copy data and launch kernels in multiple concurrent streams.
    for (int k = 0; k < howmany; k++) {
	const cudaStream_t stream = strm[k % MAX_STRM];
	memcpy(pinned + k*dist,  Qlm + (k%MAX_STRM)*nlm, 2*nlm*sizeof(double));	// copy to pinned mem
	cudaMemcpyAsync(d_qlm + (k%MAX_STRM)*2*nlm, pinned + k*dist, 2*nlm*sizeof(double), cudaMemcpyHostToDevice, stream);		// data transfer to gpu
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (mmax==0) {
	    leg_m0<0><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, d_qlm + (k%MAX_STRM)*2*nlm, d_q + (k%MAX_STRM)*nspat, llim, nlat/2);
	} else if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<0><<<blocks, threads, 2*threadsPerBlock*sizeof(double), stream>>>(d_alm, d_ct, d_qlm + (k%MAX_STRM)*2*nlm, d_q + (k%MAX_STRM)*nspat, llim, nlat/2, lmax,mres, nphi);
	} else {
	    leg_m_highllim<0><<<blocks, threads, 0, stream>>>(d_alm, d_ct, d_qlm + (k%MAX_STRM)*2*nlm, d_q + (k%MAX_STRM)*nspat, llim, nlat/2, lmax,mres, nphi);
	}
	// padd missing m's with 0 (m>mmax)
	if (2*(mmax+1) <= nphi)
	    cudaMemsetAsync( d_q + (k%MAX_STRM)*nspat + (mmax+1)*nlat, 0, sizeof(double)*(nphi-2*mmax-1)*nlat, stream );		// set to zero before fft

	cudaMemcpyAsync(pinned + k*dist, d_q + (k%MAX_STRM)*nspat, nspat*sizeof(double), cudaMemcpyDeviceToHost, stream);		// transfer back
	cudaEventCreateWithFlags(event + k,  cudaEventBlockingSync | cudaEventDisableTiming);
	cudaEventRecord(event[k], stream);
    }

    // wait for event completion, and perform fft (out-of-place) by the CPU (possibly multi-threaded).
    for (int k = 0; k < howmany; k++) {
	cudaEventSynchronize( event[k] );		// wait for event
	cudaEventDestroy(event[k]);			// get rid of event
	if (nphi > 1) {		// fft
	    if (shtns->fftc_mode == 0) {
		memcpy(Vr + (k%MAX_STRM)*nspat, pinned + k*dist, nspat*sizeof(double));	// copy from pinned mem (should be replaced by oop fft)
		fftw_execute_dft(shtns->ifftc, (cplx *) (Vr + (k%MAX_STRM)*nspat), (cplx *) (Vr + (k%MAX_STRM)*nspat));
		//fftw_execute_dft(shtns->ifftc, (cplx *) (pinned + k*dist), (cplx *) (Vr + k*nspat));	// oop
	    } else {		// split dft
		printf("ERROR fft not supported\n");
	    }
	}
    }

    free(event);
    cudaFreeHost(pinned);
}

