function [n] = fft_int(m, fmax)

if (m <= fmax)
  n=m;
  return;
end;

if (fmax <= 1)
  n=0;
  return;
end;

n=m-1;
if (mod(n,2) > 0)
  n = n-1;	% even only
end;

f=0;
while (f != n)

n = n+2;

f=1;
for k=2:fmax
  while (mod(n,k*f) == 0) & (k*f <=n)
    f = k*f;
  end;
end;

end;

