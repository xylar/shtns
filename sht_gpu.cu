
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include "sht_private.h"

#include <cufft.h>

// 256 for scalar SH_to_spat seems best on kepler.
#define THREADS_PER_BLOCK 128

// adjustment for cuda
#undef SHT_L_RESCALE_FLY
#undef SHT_ACCURACY
#define SHT_L_RESCALE_FLY 1800
#define SHT_ACCURACY 1.0e-40


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
    const double cost = (it<nlat_2) ? ct[it] : 0.0;
    double y0 = ak[0];
    double re = y0 * qk[0];
    double y1 = y0 * ak[1] * cost;
    double ro = y1 * qk[1];
    al+=2;    l+=2;	k+=2;	kq+=2;
    while(l<llim) {
	if (k+6 > THREADS_PER_BLOCK) {
	    __syncthreads();	    
	    ak[j] = al[j];
	    if ((j&1) == 0) qk[j/2] = ql[2*l+j];
	    k=0;	kq=0;
    	    __syncthreads();
	}
	y0  = ak[k+1]*(cost*y1) + ak[k]*y0;
	re += y0 * qk[kq];
	y1  = ak[k+3]*(cost*y0) + ak[k+2]*y1;
	ro += y1 * qk[kq+1];
	al+=4;	l+=2;	k+=4;	kq+=2;
    }
    if (l==llim) {
	y0  = ak[k+1]*cost*y1 + ak[k]*y0;
	re += y0 * qk[kq];
    }

    if (it < nlat_2) {
        q[it] = re+ro;
        q[nlat_2*2-1-it] = re-ro;
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

    if (im==0) {
	ni[j] = al[j];
	si[j] = ql[j];
	__syncthreads();
	int l = 0;
	const double cost = (it < nlat_2) ? ct[it] : 0.0;
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
	    y0  = ni[ka+1]*(cost*y1) + ni[ka]*y0;
	    re += y0 * si[2*kq];
	    y1  = ni[ka+3]*(cost*y0) + ni[ka+2]*y1;
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
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*l;

	ni[j] = al[j];
	si[j] = ql[2*m+j];

	// add polar optimization ?

	const double cost = (it < nlat_2) ? ct[it] : 0.0;
	double rer,ror, rei, roi, y0, y1;
	y0 = 1.0;
	l = m;
	double stx = sqrt(1.0 - cost*cost);
	do {		// sin(theta)^m
	    if (l&1) y0 *= stx;
	    stx *= stx;
	} while(l >>= 1);

	__syncthreads();
	y0 *= ni[0];
	ror = 0.0;		roi = 0.0;
	rer = 0.0;		rei = 0.0;
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

    __shared__ double ni[THREADS_PER_BLOCK];
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

	// add polar optimization

	double rer,ror, rei, roi, y0, y1;
	ni[j] = 1.0;	// y0
	l = m;
	int ny = 0;
	si[j] = sqrt(1.0 - cost*cost);	// stx
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
	ror = 0.0;		roi = 0.0;
	rer = 0.0;		rei = 0.0;
	y1 = al[1]*y0*cost;

	l=m;		al+=2;
	if (ny<0) {
	    ni[j] = al[j&31];
	    int ka = 0;
	    while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
		if (ka+4 > 32) {
		    ni[j] = al[(j&31)];
		    ka=0;
		}
		//y0 = al[1]*cost*y1 + al[0]*y0;
		//y1 = al[3]*cost*y0 + al[2]*y1;
		y0 = ni[ka+1+(j&0xFFE0)]*(cost*y1) + ni[ka+(j&0xFFE0)]*y0;
		y1 = ni[ka+3+(j&0xFFE0)]*(cost*y0) + ni[ka+2+(j&0xFFE0)]*y1;
		l+=2;	al+=4;	ka+=4;
		si[j] = y0;
		if (fabs(si[(j&0xE0)+31]) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
		    ++ny;
		    y0 *= 1.0/SHT_SCALE_FACTOR;
		    y1 *= 1.0/SHT_SCALE_FACTOR;
		}
	    }
	}
	if (ny == 0) {
	    ni[j] = al[j&31];
	    si[j] = ql[2*l+(j&31)];
	    int kq = 0;		int ka = 0;
	    while (l<llim) {	// compute even and odd parts
		if (2*kq+4 > 32) {
		    ni[j] = al[(j&31)];
		    si[j] = ql[2*l+(j&31)];
		    ka=0;	kq=0;
		}
		rer += y0 * si[2*kq+(j&0xFFE0)];	// real
		rei += y0 * si[2*kq+1+(j&0xFFE0)];	// imag
		y0 = ni[ka+1+(j&0xFFE0)]*(cost*y1) + ni[ka+(j&0xFFE0)]*y0;
		ror += y1 * si[2*kq+2+(j&0xFFE0)];	// real
		roi += y1 * si[2*kq+3+(j&0xFFE0)];	// imag
		y1 = ni[ka+3+(j&0xFFE0)]*(cost*y0) + ni[ka+2+(j&0xFFE0)]*y1;
		l+=2;	al+=4;	kq+=2;	ka+=4;
	    }
	    if (l==llim) {
		rer += y0 * ql[2*l];
		rei += y0 * ql[2*l+1];
	    }
	}

	if (it < nlat_2) {
	    /// store mangled for complex fft
	    // first, we store to shared memory the north and south values
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
int shtns_init_gpu(shtns_cfg shtns)
{
    cudaError_t err = cudaSuccess;
    const long nlm = shtns->nlm;
    const long nlat_2 = shtns->nlat_2;

    double *d_alm = NULL;
    double *d_ct  = NULL;
    double *d_qlm = NULL;
    double *d_q   = NULL;
    int err_count = 0;

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, 0);
    if (err != cudaSuccess) return 0;
    #if SHT_VERBOSE > 0
    printf("  cuda GPU \"%s\" found (warp size = %d).\n", prop.name, prop.warpSize);
    #endif
    if (prop.warpSize != 32) return 0;		// failure, SHTns requires a warpSize of 32.

    // Allocate the device input vector alm
    err = cudaMalloc((void **)&d_alm, 2*nlm*sizeof(double));
    if (err != cudaSuccess) err_count ++;
    // Allocate the device input vector ct
    err = cudaMalloc((void **)&d_ct, nlat_2*sizeof(double));
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
	return 0;	// fail
    }

    shtns->d_alm = d_alm;
    shtns->d_ct  = d_ct;
    shtns->d_q   = d_q;
    shtns->d_qlm = d_qlm;
    return 1;	// success
}

extern "C"
void shtns_release_gpu(shtns_cfg shtns)
{
    cufftDestroy(shtns->cufft_plan);
    cudaFree(shtns->d_q);
    cudaFree(shtns->d_qlm);
    cudaFree(shtns->d_ct);
    cudaFree(shtns->d_alm);
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
    const int blocksPerGrid =(nlat/2 + threadsPerBlock - 1) / threadsPerBlock;
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
    if (err != cudaSuccess) { printf("CUDA error : %s!\n", cudaGetErrorString(err));	return; }

    // copy back spatial data
    err = cudaMemcpy(Vr, d_q, nlat*nphi*sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { printf("failed copy back\n");	return; }
}

