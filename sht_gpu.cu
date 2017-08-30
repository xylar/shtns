
/* TODO
 * 1) try to store data for complex2real fft, and perform fft on host (less data to transfer)
 * 2) use static polar optimization (from constant memory ?)
 * 3) use a for loop in m-direction to re-use threads at larger m's.
 * 4) use double2 (double4 ?) vector types to avoid the use of shared memory to write mangled data.
 * 5) multi-stream / multi-gpu computation ?
 */

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include "sht_private.h"

#include <cufft.h>

// 256 for scalar SH_to_spat seems best on kepler.
#define THREADS_PER_BLOCK 256
// number of latitudes per thread:
#define NWAY 2
// the warp size is always 32 on cuda devices (up to Pascal at least)
#define WARPSZE 32

// adjustment for cuda
#undef SHT_L_RESCALE_FLY
#undef SHT_ACCURACY
#define SHT_L_RESCALE_FLY 1800
#define SHT_ACCURACY 1.0e-40


/// On KEPLER, This kernel is fastest with THREADS_PER_BLOCK=256 and NWAY=1
__global__ void
leg_m0(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2)
{
    // im = 0
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int j = threadIdx.x;

    __shared__ double ak[THREADS_PER_BLOCK];
    __shared__ double qk[THREADS_PER_BLOCK/2];
    ak[j] = al[j];
    if ((j&1) == 0) qk[j/2] = ql[j];
    __syncthreads();

    int l = 0;
    int k = 0;	int kq = 0;
    double cost[NWAY];
    double y0[NWAY];    double y1[NWAY];
    double re[NWAY];    double ro[NWAY];

    for (int i=0; i<NWAY; i++) {
	cost[i] = (it+i<nlat_2) ? ct[it+i] : 0.0;
	y0[i] = ak[0];
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
	if (k+6 > THREADS_PER_BLOCK) {
	    __syncthreads();
	    ak[j] = al[j];
	    if ((j&1) == 0) qk[j/2] = ql[2*l+j];
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
void warp_reduce_add_2(double& ev, double& od) {
  for (int offset = WARPSZE/2; offset > 0; offset >>= 1) {
    ev += __shfl_down(ev, offset);
    od += __shfl_down(od, offset);
  }
}

__inline__ __device__
void warp_reduce_add(double& ev) {
  for (int offset = WARPSZE/2; offset > 0; offset >>= 1) {
    ev += __shfl_down(ev, offset);
  }
}

#if __CUDA_ARCH__ < 600
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

/// On KEPLER, works best with THREADS_PER_BLOCK=32 and NWAY=3
__global__ void
ileg_m0(const double *al, const double *ct, const double *q, double *ql, const int llim, const int nlat_2)
{
    // im = 0
    const int it = (blockDim.x * blockIdx.x + threadIdx.x)*NWAY;
    const int j = threadIdx.x;

    __shared__ double ak[THREADS_PER_BLOCK];
    #if THREADS_PER_BLOCK > WARPSZE
    __shared__ double qk[2*THREADS_PER_BLOCK/WARPSZE];
    #endif

    double y0[NWAY], y1[NWAY];
    double re[NWAY],ro[NWAY],  qe,qo;
    for (int i=0; i<NWAY; i++) {
	y0[i] = (it+i < nlat_2) ? q[it+i] : 0.0;		// north
	y1[i] = (it+i < nlat_2) ? q[nlat_2*2-1-(it+i)] : 0.0;	// south
	re[i] = y0[i] + y1[i];
	ro[i] = y0[i] - y1[i];
    }

    ak[j] = al[j];
    #if THREADS_PER_BLOCK > WARPSZE
    __syncthreads();
    #endif

    int l = 0;
    double cost[NWAY];
    for (int i=0; i<NWAY; i++) {
	if (it+i < nlat_2) {
	    y0[i] = ct[it+i + nlat_2];		// weight are stored just after ct.
	    cost[i] = ct[it+i];
	} else {
	    y0[i] = 0.0;	    cost[i] = 0.0;
	}
    }
    for (int i=0; i<NWAY; i++) 	y0[i] *= ak[0];
    for (int i=0; i<NWAY; i++) 	y1[i] = y0[i] * ak[1] * cost[i];
    al+=2;	int k=2;
    while (l<llim) {
	qe = re[0]*y0[0];
	qo = ro[0]*y1[0];
	y0[0]  = ak[k+1]*cost[0]*y1[0] + ak[k]*y0[0];
	y1[0]  = ak[k+3]*cost[0]*y0[0] + ak[k+2]*y1[0];
	for (int i=1; i<NWAY; i++) {
	    qe += re[i]*y0[i];
	    qo += ro[i]*y1[i];
	    y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
	    y1[i]  = ak[k+3]*cost[i]*y0[i] + ak[k+2]*y1[i];
	}
	al+=4;		k+=4;	// 4 als consumed.
	warp_reduce_add_2(qe, qo);		// reduce within warp
	#if 2*THREADS_PER_BLOCK > WARPSZE*WARPSZE
	    #error "threads/block should be <= warp-size^2"
	#endif
	#if THREADS_PER_BLOCK > WARPSZE
	if ((j&31) == 0) {			// first lane stores in shared mem
	    qk[2*j/WARPSZE] = qe;
	    qk[2*j/WARPSZE+1] = qo;
	}
	__syncthreads();
	qe = (j<2*THREADS_PER_BLOCK/WARPSZE) ? qk[j] : 0.0;			// read from shared mem
	for (int offset = THREADS_PER_BLOCK/WARPSZE; offset > 1; offset >>= 1)
	    qe += __shfl_down(qe, offset);
	#else
	qo = __shfl_up(qo, 1,2);
	if (j==1) qe = qo;
	#endif
	if (k+4 > THREADS_PER_BLOCK) {		// cache next als (between the two __syncthreads() calls)
	    ak[j] = al[j];	k=0;
	}
	#if THREADS_PER_BLOCK > WARPSZE
	__syncthreads();
	#endif
	if (j<2) {		// should be an atomic-add to global mem for large transforms (nlat_2 > block_size)
	    if (nlat_2 <= THREADS_PER_BLOCK) {		// do we need atomic add or not ?
		ql[2*(l+j)] = qe;
//		ql[2*(l+j)+1] = 0.0;
	    } else {
		atomicAdd(ql+2*(l+j), qe);		// VERY slow atomic add on Kepler.
	    }
	}
	l+=2;
    }
    if (l==llim) {
	qe = re[0]*y0[0];
	for (int i=1; i<NWAY; i++) {
	    qe += re[i]*y0[i];
	}
	warp_reduce_add(qe);
	#if THREADS_PER_BLOCK > WARPSZE
	if ((j&31) == 0) {			// first lane stores in shared mem
	    qk[j/WARPSZE] = qe;
	}
	__syncthreads();	// block-reduce, max block-size = 1024
	qe = (j<THREADS_PER_BLOCK/WARPSZE) ? qk[j] : 0.0;			// read from shared mem
	for (int offset = THREADS_PER_BLOCK/WARPSZE/2; offset > 0; offset >>= 1)
	    qe += __shfl_down(qe, offset);
	#endif
	if (j==0) {		// should be an atomic-add to global mem for large transforms (nlat_2 > block_size)
	    if (nlat_2 <= THREADS_PER_BLOCK) {
	        ql[2*l] = qe;
//	        ql[2*l+1] = 0.0;
	    } else {
		atomicAdd(ql+2*l, qe);		// VERY slow atomic add on Kepler.
	    }
	}
    }
}



/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim MUST BE <= 1800
__global__ void
leg_m_lowllim(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2, const int lmax, const int mmax, const int mres, const int nphi)
{
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int im = blockIdx.y;
    const int j = threadIdx.x;
    const int m_inc = 2*nlat_2;
    const int k_inc = 1;

    __shared__ double ni[THREADS_PER_BLOCK];
    __shared__ double si[THREADS_PER_BLOCK];

    const double cost = (it < nlat_2) ? ct[it] : 0.0;

    if (im==0) {
	ni[j] = al[j];
	si[j] = ql[j];
	__syncthreads();
	int l = 0;
	int ka = 0;	int kq = 0;
	double y0 = ni[0];
	double re = y0 * si[0];
	double y1 = y0 * ni[1] * cost;
	double ro = y1 * si[2];
	al+=2;    l+=2;		ka+=2;	kq+=2;
	while(l<llim) {
	    if (ka+6 > THREADS_PER_BLOCK) {
		__syncthreads();  
		ni[j] = al[j];
		si[j] = ql[2*l+j];
		ka=0;	kq=0;
		__syncthreads();
	    }
	    y0  = ni[ka+1]*cost*y1 + ni[ka]*y0;
	    re += y0 * si[2*kq];
	    y1  = ni[ka+3]*cost*y0 + ni[ka+2]*y1;
	    ro += y1 * si[2*kq+2];
	    al+=4;	l+=2;	  ka+=4;    kq+=2;
	}
	if (l==llim) {
	    y0  = ni[ka+1]*cost*y1 + ni[ka]*y0;
	    re += y0 * si[2*kq];
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
	ql += 2*l;

	y1 = sqrt(1.0 - cost*cost);		// y1 = sin(theta)
	ni[j] = al[j];
	si[j] = ql[2*m+j];

	ror = 0.0;		roi = 0.0;
	rer = 0.0;		rei = 0.0;
	y0 = 1.0;
	l = m;
	do {		// sin(theta)^m
	    if (l&1) y0 *= y1;
	    y1 *= y1;
	} while(l >>= 1);

	__syncthreads();
	y0 *= ni[0];
	y1 = ni[1]*y0*cost;

	int ka = 2;
	l=m;		al+=2;
	int kq = 0;

	while (l<llim) {	// compute even and odd parts
	    if (2*kq+6 > THREADS_PER_BLOCK) {
		__syncthreads();
		ni[j] = al[j];
		si[j] = ql[2*l+j];
		ka=0;	kq=0;
		__syncthreads();
	    }
	    rer += y0 * si[2*kq];	// real
	    rei += y0 * si[2*kq+1];	// imag
	    y0 = ni[ka+1]*(cost*y1) + ni[ka]*y0;
	    ror += y1 * si[2*kq+2];	// real
	    roi += y1 * si[2*kq+3];	// imag
	    y1 = ni[ka+3]*(cost*y0) + ni[ka+2]*y1;
	    l+=2;	al+=4;	 ka+=4;	  kq+=2;
	}
	if (l==llim) {
	    rer += y0 * si[2*kq];
	    rei += y0 * si[2*kq+1];
	}

	__syncthreads();
	if (it < nlat_2) {
	    /// store mangled for complex fft
	    // first, we store to shared memory the north and south values
	    double nr = rer+ror;
	    double sr = rer-ror;
	    #if __CUDA_ARCH__ < 300
	    ni[j] = rei+roi;
	    si[j] = rei-roi;
	    // __syncthreads();
	    // combine and store to global memory.
	    const int xchg = 1 - 2*(j&1);		// +/- 1
	    const double sgn = xchg;
	    double nix = sgn * ni[j+xchg];
	    double six = sgn * si[j+xchg];
	    #else
	    const double sgn = 1 - 2*(j&1);
	    rei = __shfl_xor(rei, 1);
	    roi = __shfl_xor(roi, 1);
	    double nix = sgn*(rei+roi);
	    double six = sgn*(rei-roi);
	    #endif	    
	    q[im*m_inc + it*k_inc]                     = nr - nix;
	    q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
	    q[im*m_inc + (nlat_2*2-1-it)*k_inc]        = sr + six;
	    q[(nphi-im)*m_inc + (nlat_2*2-1-it)*k_inc] = sr - six;

	    if (im==mmax) {		// write extra zeros (padding for m>mmax)
		int iim = im+1;
		while(iim <= nphi-(mmax+1)) {	// padding for high m's
			q[iim*m_inc + it*k_inc] = 0.0;
			q[iim*m_inc + (nlat_2*2-1-it)*k_inc] = 0.0;
			iim += 1;
		}
	    }
	}
    }
}

/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim can be arbitrarily large (> 1800)
__global__ void
leg_m_highllim(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2, const int lmax, const int mmax, const int mres, const int nphi)
{
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int im = blockIdx.y;
    const int j = threadIdx.x;
    const int m_inc = 2*nlat_2;
    const int k_inc = 1;

    __shared__ double ni[THREADS_PER_BLOCK];	// cache
    __shared__ double si[THREADS_PER_BLOCK];

    const double cost = (it < nlat_2) ? ct[it] : 0.0;

    if (im==0) {
	int l = 0;
	double y0 = al[0];
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
	ql += 2*l;
	double rer,ror, rei,roi, y0, y1;
	ror = 0.0;	roi = 0.0;
	rer = 0.0;	rei = 0.0;
	#if __CUDA_ARCH__ < 300
	si[j] = sqrt(1.0 - cost*cost);	// sin(theta)
	if (m - llim*si[(j&0xE0) +31] <= max(50,llim/200)) {		// polar optimization (see Reinecke 2013), avoiding warp divergence
	    ni[j] = 1.0;	// y0
	    l = m;
	    int ny = 0;
	    int nsint = 0;
	    do {		// sin(theta)^m		(use rescaling to avoid underflow)
		if (l&1) {
		    ni[j] *= si[j];
		    ny += nsint;
		    if (ni[(j&0xE0) +31] < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {		// avoid warp divergence here !
			ny--;
			ni[j] *= SHT_SCALE_FACTOR;
		    }
		}
		si[j] *= si[j];
		nsint += nsint;
		if (si[(j&0xE0) +31] < 1.0/SHT_SCALE_FACTOR) {		// avoid warp divergence here !
		    nsint--;
		    si[j] *= SHT_SCALE_FACTOR;
		}
	    } while(l >>= 1);
	    y0 = ni[j] * al[0];
	    y1 = al[1]*y0*cost;
	#else
	y1 = sqrt(1.0 - cost*cost);	// sin(theta)
	if (m - llim*__shfl(y1,31) <= max(50,llim/200)) {		// polar optimization (see Reinecke 2013), avoiding warp divergence
	    y0 = 1.0;	// y0
	    l = m;
	    int ny = 0;
	    int nsint = 0;
	    do {		// sin(theta)^m		(use rescaling to avoid underflow)
		if (l&1) {
		    y0 *= y1;
		    ny += nsint;
		    if (__shfl(y0, 31) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {		// avoid warp divergence
			ny--;
			y0 *= SHT_SCALE_FACTOR;
		    }
		}
		y1 *= y1;
		nsint += nsint;
		if (__shfl(y1, 31) < 1.0/SHT_SCALE_FACTOR) {		// avoid warp divergence
		    nsint--;
		    y1 *= SHT_SCALE_FACTOR;
		}
	    } while(l >>= 1);
	    y0 *= al[0];
	    y1 = al[1]*y0*cost;
	#endif

	    l=m;		al+=2;
	    if (ny<0) {
		ni[j] = al[j&31];
		int ka = 0;
		while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
		    if (ka+4 > WARPSZE) {
			ni[j] = al[(j&31)];
			ka=0;
		    }
		    //y0 = al[1]*cost*y1 + al[0]*y0;
		    //y1 = al[3]*cost*y0 + al[2]*y1;
		    y0 = ni[ka+1+(j&0xFFE0)]*(cost*y1) + ni[ka+(j&0xFFE0)]*y0;
		    y1 = ni[ka+3+(j&0xFFE0)]*(cost*y0) + ni[ka+2+(j&0xFFE0)]*y1;
		    l+=2;	al+=4;	ka+=4;
		    #if __CUDA_ARCH__ < 300
		    si[j] = y0;
		    if (fabs(si[(j&0xE0)+31]) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0)
		    #else
		    if (fabs(__shfl(y0,31)) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0)
		    #endif
		    {	// rescale when value is significant
			++ny;
			y0 *= 1.0/SHT_SCALE_FACTOR;
			y1 *= 1.0/SHT_SCALE_FACTOR;
		    }
		}
	    }
	    if (ny == 0) {
		ni[j] = al[j&31];
		si[j] = ql[2*l+(j&31)];
		int k = 0;		const int ofs = j & 0xFFE0;
		while (l<llim) {	// compute even and odd parts
		    if (k+4 > WARPSZE) {
			ni[j] = al[(j&31)];
			si[j] = ql[2*l+(j&31)];
			k=0;
		    }
		    rer += y0 * si[k+ofs];	// real
		    rei += y0 * si[k+1+ofs];	// imag
		    y0 = ni[k+1+ofs]*(cost*y1) + ni[k+ofs]*y0;
		    ror += y1 * si[k+2+ofs];	// real
		    roi += y1 * si[k+3+ofs];	// imag
		    y1 = ni[k+3+ofs]*(cost*y0) + ni[k+2+ofs]*y1;
		    l+=2;	al+=4;	k+=4;
		}
		if (l==llim) {
		    rer += y0 * ql[2*l];
		    rei += y0 * ql[2*l+1];
		}
	    }
	}

	if (it < nlat_2) {
	    /// store mangled for complex fft
	    double nr = rer+ror;
	    double sr = rer-ror;
	    #if __CUDA_ARCH__ < 300
	    ni[j] = rei+roi;
	    si[j] = rei-roi;
	    const int xchg = 1 - 2*(j&1);		// +/- 1
	    const double sgn = xchg;
	    double nix = sgn * ni[j+xchg];
	    double six = sgn * si[j+xchg];	    
	    #else
	    const double sgn = 1 - 2*(j&1);
	    rei = __shfl_xor(rei, 1);
	    roi = __shfl_xor(roi, 1);
	    double nix = sgn*(rei+roi);
	    double six = sgn*(rei-roi);
	    #endif
	    q[im*m_inc + it*k_inc]                     = nr - nix;
	    q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
	    q[im*m_inc + (nlat_2*2-1-it)*k_inc]        = sr + six;
	    q[(nphi-im)*m_inc + (nlat_2*2-1-it)*k_inc] = sr - six;

	    if (im==mmax) {		// write extra zeros (padding for m>mmax)
		int iim = im+1;
		while(iim <= nphi-(mmax+1)) {	// padding for high m's
			q[iim*m_inc + it*k_inc] = 0.0;
			q[iim*m_inc + (nlat_2*2-1-it)*k_inc] = 0.0;
			iim += 1;
		}
	    }
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
    int err_count = 0;
    int device_id = -1;

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, 0);
    if (err != cudaSuccess) return -1;
    #if SHT_VERBOSE > 0
    printf("  cuda GPU \"%s\" found (warp size = %d, compute capabilities = %d.%d).\n", prop.name, prop.warpSize, prop.major, prop.minor);
    #endif
    if (prop.warpSize != WARPSZE) return -1;		// failure, SHTns requires a warpSize of 32.

    // Allocate the device input vector alm
    err = cudaMalloc((void **)&d_alm, 2*nlm*sizeof(double));
    if (err != cudaSuccess) err_count ++;
    // Allocate the device input vector cos(theta) and gauss weights
    err = cudaMalloc((void **)&d_ct, 2*nlat_2*sizeof(double));
    if (err != cudaSuccess) err_count ++;
    // Allocate the device work vector qlm
    err = cudaMalloc((void **)&d_qlm, 2*nlm*sizeof(double));
    if (err != cudaSuccess) err_count ++;
    // Allocate the device work vector q
    err = cudaMalloc((void **)&d_q, shtns->nlat * shtns->nphi * sizeof(double));
    if (err != cudaSuccess) err_count ++;

    if (err_count == 0) {
	err = cudaMemcpy(d_alm, shtns->alm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
	err = cudaMemcpy(d_ct, shtns->ct, nlat_2*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
	err = cudaMemcpy(d_ct + nlat_2, shtns->wg, nlat_2*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  err_count ++;
    }

    /* cuFFT init */
    int nfft = shtns->nphi;
    if (nfft > 1) {
	// cufftPlanMany(cufftHandle *plan, int rank, int *n,   int *inembed, int istride, int idist,   int *onembed, int ostride, int odist,   cufftType type, int batch);
	cufftResult res;
	res = cufftPlanMany((cufftHandle*) &shtns->cufft_plan, 1, &nfft, &nfft, shtns->nlat_2, 1, &nfft, shtns->nlat_2, 1, CUFFT_Z2Z, shtns->nlat_2);
	if (res != CUFFT_SUCCESS)  err_count ++;
    }

    if (err_count != 0) {
	cudaFree(d_q);	cudaFree(d_qlm);  cudaFree(d_ct);  cudaFree(d_alm);
	return -1;	// fail
    }

    shtns->d_alm = d_alm;
    shtns->d_ct  = d_ct;
    shtns->d_q   = d_q;
    shtns->d_qlm = d_qlm;
    cudaGetDevice(&device_id);
    return device_id;		// success, return device_id
}

extern "C"
void cushtns_release_gpu(shtns_cfg shtns)
{
    cufftDestroy(shtns->cufft_plan);
    cudaFree(shtns->d_q);
    cudaFree(shtns->d_qlm);
    cudaFree(shtns->d_ct);
    cudaFree(shtns->d_alm);
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
    const int mmax = shtns->mmax;
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
    if (mmax == 0) {
	double* Ql0;
	Ql0 = (double*) malloc((lmax+1)*sizeof(double));
	for (int l=0; l<=llim; l++) {
	    Ql0[l] = creal(Qlm[l]);
	}
	err = cudaMemcpy(d_qlm, Ql0, (llim+1)*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	leg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2);
    } else {
	err = cudaMemcpy(d_qlm, Qlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<<<blocks, threads>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2, lmax,mmax,mres, nphi);
	} else {
	    leg_m_highllim<<<blocks, threads>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2, lmax,mmax,mres, nphi);
	}
    }
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Failed to launch cuda kernel (error code %s)!\n", cudaGetErrorString(err));
    }

    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)  printf("failed copy back\n");

    if (mmax > 0) {		// fft
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
extern "C"
void cu_SH_to_spat(shtns_cfg shtns, cplx* d_Qlm, double *d_Vr, const long int llim)
{
    const int lmax = shtns->lmax;
    const int mmax = shtns->mmax;
    const int mres = shtns->mres;
    const int nlat_2 = shtns->nlat_2;
    const int nphi = shtns->nphi;
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat_2 + threadsPerBlock - 1) / threadsPerBlock;
    if (nphi == 1) {
	leg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2);
    } else {
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (mmax==0) {
	    leg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2);
	} else
	if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<<<blocks, threads>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mmax,mres, nphi);
	} else {
	    leg_m_highllim<<<blocks, threads>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mmax,mres, nphi);
	}
	cufftResult res;
	res = cufftExecZ2Z((cufftHandle) shtns->cufft_plan, (cufftDoubleComplex*) d_Vr, (cufftDoubleComplex*) d_Vr, CUFFT_INVERSE);
	if (res != CUFFT_SUCCESS) printf("cufft error %d\n", res);
    }
}

/// Perform SH transform on data that is already on the GPU. d_Qlm and d_Vr are pointers to GPU memory (obtained by cudaMalloc() for instance)
extern "C"
void cu_spat_to_SH(shtns_cfg shtns, double *d_Vr, cplx* d_Qlm, const long int llim)
{
    const int mmax = shtns->mmax;
    const int mres = shtns->mres;
    const int nlat_2 = shtns->nlat_2;
    const int nphi = shtns->nphi;
    const int nlm = shtns->nlm;
    double *d_alm = shtns->d_alm;
    double *d_ct = shtns->d_ct;

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat_2 + threadsPerBlock*NWAY - 1) / (threadsPerBlock*NWAY);
    cudaMemset(d_Qlm, 0, sizeof(double)*2*nlm);		// set to zero before we start.
    if (nphi == 1) {
	ileg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2);
    } else {
	cufftResult res;
	res = cufftExecZ2Z((cufftHandle) shtns->cufft_plan, (cufftDoubleComplex*) d_Vr, (cufftDoubleComplex*) d_Vr, CUFFT_INVERSE);
	if (res != CUFFT_SUCCESS) printf("cufft error %d\n", res);

	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	if (mmax==0) {
	    ileg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, (double*) d_Vr, (double*) d_Qlm, llim, nlat_2);
	} else {
	    printf("NOT implemented\n");	return;
	}
/*	if (llim <= SHT_L_RESCALE_FLY) {
	    leg_m_lowllim<<<blocks, threads>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mmax,mres, nphi);
	} else {
	    leg_m_highllim<<<blocks, threads>>>(d_alm, d_ct, (double*) d_Qlm, (double*) d_Vr, llim, nlat_2, lmax,mmax,mres, nphi);
	}
	*/
    }
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

    // copy spectral data to GPU
    err = cudaMemcpy(d_qlm, Qlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { printf("failed copy qlm\n");	return; }

    // SHT on the GPU
    cu_SH_to_spat(shtns, (cplx*) d_qlm, d_q, llim);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("SH_to_spat CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spatial data
    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { printf("failed copy back\n");	return; }
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
    if (err != cudaSuccess) { printf("failed copy q\n");	return; }

    // SHT on the GPU
    cu_spat_to_SH(shtns, d_q, (cplx*) d_qlm, llim);
    err = cudaGetLastError();
    if (err != cudaSuccess) { printf("spat_to_SH CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spectral data
    err = cudaMemcpy(Qlm, d_qlm, 2*nlm*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { printf("failed copy back\n");	return; }
}

