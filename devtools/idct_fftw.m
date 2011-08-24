% dct qui se comporte comme celle de fftw

function z = idct_fftw( dz )

sz=size(dz);
N = sz(1);
if (N == 1) 
  N = sz(2);
  dz(:,1) = dz(:,1) * 2*sqrt(2);
  z = idct(dz) * 0.5*sqrt(N/2);
else
  dz(1,:) = dz(1,:) * 2*sqrt(2);
  z = idct(dz) * 0.5*sqrt(N/2);
end

