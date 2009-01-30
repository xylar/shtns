% dct qui se comporte comme celle de fftw

function dz = dct_fftw( z )

sz=size(z);
N = sz(1);
if (N == 1) 
  N = sz(2);
  dz = dct(z) * sqrt(2*N);
  dz(:,1) = dz(:,1) * sqrt(2);
else
  dz = dct(z) * sqrt(2*N);
  dz(1,:) = dz(1,:) * sqrt(2);
end

