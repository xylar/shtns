

i=1;
n=1;

while (n<4096)
  n = fft_int(n+1,7);
  ilist(i) = n;
  i = i+1;
end;
