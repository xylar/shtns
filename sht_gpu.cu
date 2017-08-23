
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include "sht_private.h"

#define THREADS_PER_BLOCK 32


__global__ void
leg_m0(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2)
{
    // im = 0
    const int it = blockDim.x * blockIdx.x + threadIdx.x;

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
}

__global__ void
leg_m(const double *al, const double *ct, const double *ql, double *q, const int llim, const int nlat_2, const int lmax, const int mmax, const int mres, const int nphi)
{
    const int it = blockDim.x * blockIdx.x + threadIdx.x;
    const int im = blockIdx.y;
    const int m_inc = 2*nlat_2;
    const int k_inc = 1;

    if (im==0) {
	if (it < nlat_2) {
	    int l = 0;
	    const double cost = ct[it];
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
	    // store mangled for complex fft
	    q[it*k_inc] = re+ro;
	    q[(nlat_2*2-1-it)*k_inc] = re-ro;
	}
    } else { 	// m>0
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*l;

	// add polar optimization
	if (it < nlat_2) {
	    const double cost = ct[it];
	    double rer,ror, rei, roi, y0, y1;
	    y0 = 1.0;
	    int l = im*mres;
	    int ny = 0;
	    double stx = sqrt(1.0 - cost*cost);
	    if (llim <= SHT_L_RESCALE_FLY) {
		do {		// sin(theta)^m
		    if (l&1) y0 *= stx;
		    stx *= stx;
		} while(l >>= 1);
	    } else {
		long int nsint = 0;
		do {		// sin(theta)^m		(use rescaling to avoid underflow)
		    if (l&1) {
			y0 *= stx;
			ny += nsint;
			if (y0 < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {		// TODO avoid thread divergence here !
			    ny--;
			    y0 *= SHT_SCALE_FACTOR;
			}
		    }
		    stx *= stx;
		    nsint += nsint;
		    if (stx < 1.0/SHT_SCALE_FACTOR) {
			nsint--;
			stx *= SHT_SCALE_FACTOR;
		    }
		} while(l >>= 1);
	    }
	    
	    y0 *= al[0];
	    ror = 0.0;		roi = 0.0;
	    rer = 0.0;		rei = 0.0;
	    y1 = al[1]*y0*cost;

	    l=im*mres;		al+=2;
	    while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
		y0 = al[1]*cost*y1 + al[0]*y0;
		y1 = al[3]*cost*y0 + al[2]*y1;
		l+=2;	al+=4;
		if (y0 > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
		    ++ny;
		    y0 *= 1.0/SHT_SCALE_FACTOR;
		    y1 *= 1.0/SHT_SCALE_FACTOR;
		}
	    }
	    if (ny == 0) {
		while (l<llim) {	// compute even and odd parts
		    rer += y0 * ql[2*l];	// real
		    rei += y0 * ql[2*l+1];	// imag
		    y0 = al[1]*(cost*y1) + al[0]*y0;
		    ror += y1 * ql[2*l+2];	// real
		    roi += y1 * ql[2*l+3];	// imag
		    y1 = al[3]*(cost*y0) + al[2]*y1;
		    l+=2;	al+=4;
		}
		if (l==llim) {
		    rer += y0 * ql[2*l];
		    rei += y0 * ql[2*l+1];
		}
	    }

	    /// store mangled for complex fft
	    // first, we store to shared memory the north and south values
	    __shared__ double ni[THREADS_PER_BLOCK];
	    __shared__ double si[THREADS_PER_BLOCK];
	    const int j = threadIdx.x;
	    double nr = rer+ror;
	    double sr = rer-ror;
	    ni[j] = rei+roi;
	    si[j] = rei-roi;

	    // __syncthreads();
	    // combine and store to global memory.
	    const int xchg = 1 - 2*(j&1);		// +/- 1
	    const double sgn = xchg;
	    double nix = sgn * ni[j+xchg];
	    q[im*m_inc + it*k_inc]                     = nr - nix;
	    q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
	    double six = sgn * si[j+xchg];
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
void shtns_use_gpu(shtns_cfg shtns)
{
    cudaError_t err = cudaSuccess;
    const long nlm = shtns->nlm;
    const long nlat_2 = shtns->nlat_2;

    // Allocate the device input vector alm
    double *d_alm = NULL;
    err = cudaMalloc((void **)&d_alm, 2*nlm*sizeof(double));
    if (err != cudaSuccess) printf("CUDA: failed alloc alm\n");
    // Allocate the device input vector ct
    double *d_ct = NULL;
    err = cudaMalloc((void **)&d_ct, nlat_2*sizeof(double));
    if (err != cudaSuccess) printf("CUDA: failed alloc ct\n");

    err = cudaMemcpy(d_alm, shtns->alm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)  printf("failed copy alm\n");
    err = cudaMemcpy(d_ct, shtns->ct, nlat_2*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)  printf("failed copy ct\n");    

    shtns->d_alm = d_alm;
    shtns->d_ct = d_ct;

    printf("shtns: GPU initialized\n");
}

extern "C"
void SH_to_spat_gpu(shtns_cfg shtns, cplx *Qlm, double *Vr, const long int llim)
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

    // Allocate the device output vector q
    double *d_q = NULL;
    err = cudaMalloc((void **)&d_q, nlat*nphi*sizeof(double));
    if (err != cudaSuccess) printf("failed alloc q\n");

    // Launch the Legendre CUDA Kernel
    const int threadsPerBlock = THREADS_PER_BLOCK;	// can be from 32 to 1024, we should try to measure the fastest !
    const int blocksPerGrid =(nlat/2 + threadsPerBlock - 1) / threadsPerBlock;
    double *d_qlm = NULL;
    if (mmax == 0) {
	// Allocate the device input vector qlm
	err = cudaMalloc((void **)&d_qlm, (lmax+1)*sizeof(double));
	if (err != cudaSuccess) printf("failed alloc qlm\n");
	double* Ql0;
	Ql0 = (double*) malloc((lmax+1)*sizeof(double));
	for (int l=0; l<=llim; l++) {
	    Ql0[l] = creal(Qlm[l]);
	}
	err = cudaMemcpy(d_qlm, Ql0, (llim+1)*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	leg_m0<<<blocksPerGrid, threadsPerBlock>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2);
    } else {
    
	// Allocate the device input vector qlm
	err = cudaMalloc((void **)&d_qlm, 2*nlm*sizeof(double));
	if (err != cudaSuccess) printf("failed alloc qlm\n");
	err = cudaMemcpy(d_qlm, Qlm, 2*nlm*sizeof(double), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)  printf("failed copy qlm\n");

	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	leg_m<<<blocks, threads>>>(d_alm, d_ct, d_qlm, d_q, llim, nlat/2, lmax,mmax,mres, nphi);
    }
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Failed to launch cuda kernel (error code %s)!\n", cudaGetErrorString(err));
    }

    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)  printf("failed copy back\n");

    cudaFree(d_q);
    cudaFree(d_qlm);

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
