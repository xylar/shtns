* A Fortran example program that performs backward and forward Spherical Harmonic Transforms using SHTns
      PROGRAM SHT_example

      IMPLICIT NONE
c import useful parameters for shtns initialization
      include 'shtns.f'

      Integer*4 lmax, mmax, mres
      integer*4 nlat, nphi
      integer*4 nlm
      Integer*4 layout
      complex*16, allocatable :: Slm(:), Tlm(:)
      real*8, allocatable :: Sh(:,:), Th(:,:)

      integer i,lm, m

c set size of transform
      lmax = 5
      mmax = 2
      mres = 2
      nphi = 8
      nlat = 6

c compute sizes required for arrays.
      call shtns_calc_nlm(nlm, lmax, mmax, mres)
      print*,'NLM=',nlm

c init SHT. SHT_PHI_CONTIGUOUS is defined in 'shtns.f'
      layout = SHT_PHI_CONTIGUOUS
      call shtns_init_sh_gauss(layout, lmax, mmax, mres, nlat, nphi)

c allocate memory for spectral and spatial representation
      allocate ( Slm(1:nlm), Tlm(1:nlm) )
      allocate ( Sh(1:nphi,1:nlat), Th(1:nphi,1:nlat) )
      Slm(:) = 0.0

c get index for given l and m
      call shtns_lmidx(lm, 1, 0)
      Slm(lm) = 1.0
      print*
      print*,Slm(:)

      call shtns_sh_to_spat(Slm, Sh)

c print all spatial data
      print*
      do i=1,nphi
        print*,'phi index=',i
        print*,Sh(i,:)
      enddo

      call shtns_spat_to_sh(Sh, Slm)

      print*
      print*,Slm(:)

c print SH data, grouping by m  (shows how to access SH coefficients)
      print*
      do m=0, mres*mmax, mres
        print*,'m=', m
        call shtns_lmidx(lm,m,m)
        print*,Slm(lm : lm+(lmax-m))
      enddo

      stop
      END
