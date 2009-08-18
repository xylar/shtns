* A Fortran example program that performs backward and forward Spherical Harmonic Transforms using SHTns
      PROGRAM SHT_example

      IMPLICIT NONE

      Integer*4 lmax, mmax, mres
      integer*4 nlat, nphi
      integer*4 nlm, nspat
      complex*16, allocatable :: Slm(:), Tlm(:)
      real*8, allocatable :: Sh(:), Th(:)

      integer i

c set size of transform
      lmax = 5
      mmax = 3
      mres = 1
      nphi = 8
      nlat = 6

c compute sizes required for arrays.
      call shtns_get_nlm(nlm, lmax, mmax, mres)
      print*,'NLM=',nlm
      call shtns_get_nspat_alloc(nspat, nlat, nphi)
      print*,'NSPAT_ALLOC=',nspat

c init SHT
      call shtns_init_sh_gauss(lmax, mmax, mres, nlat, nphi)

c allocate memory for spectral and spatial representation
      allocate ( Slm(1:nlm), Tlm(1:nlm) )
      allocate ( Sh(1:nspat), Th(1:nspat) )

      Slm(:) = 0.0
      Slm(2) = 1.0
      print*
      print*,Slm(:)

      call shtns_sh_to_spat(Slm, Sh)

c print all spatial data
      print*
      do i=1,nphi
        print*,'phi index=',i
        print*,Sh((i-1)*nlat+1 : i*nlat)
      enddo

      call shtns_spat_to_sh(Sh, Slm)

      print*
      print*,Slm(:)

      stop
      END
