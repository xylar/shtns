
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/// \name GPU transforms (working on GPU memory, without transfers). Non-blocking.
//@{
void cu_spat_to_SH(shtns_cfg shtns, double *Vr, cplx *Qlm, int ltr);
void cu_SH_to_spat(shtns_cfg shtns, cplx *Qlm, double *Vr, int ltr);
void cu_spat_to_SHsphtor(shtns_cfg, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int ltr);
void cu_SHsphtor_to_spat(shtns_cfg, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int ltr);
void cu_spat_to_SHqst(shtns_cfg, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr);
void cu_SHqst_to_spat(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr);
//@}

/// Set user-specified streams for compute (including fft) and transfer.
void cushtns_set_streams(shtns_cfg shtns, cudaStream_t compute_stream, cudaStream_t transfer_stream);

/// Clone a gpu-enabled shtns config, and assign it to different streams (to allow overlap and/or usage from multiple threads)
shtns_cfg cushtns_clone(shtns_cfg shtns, cudaStream_t compute_stream, cudaStream_t transfer_stream);

#ifdef __cplusplus
}
#endif /* __cplusplus */
