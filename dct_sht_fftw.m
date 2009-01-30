% test de l'algorithme...


N=10;	% point numbers
m=1;	% azimutal m
lmax=6;
ltest=4;
iter=10000;

% mapping regulier.
%t = pi*((1:N)-0.25)/(N+0.5) % nearly gauss (does NOT work)
t = pi*((1:N)-0.5)/N;	% chebychev
t = t';

% fonction test :
x = cos(t);
f = cos(t).^ltest - cos(t).^2;
df = dct_fftw(f);
dct_error = max(abs(df(ltest+2:end)))		% residual should be small

% calcule les points de Gauss
[xg,wg] = lgwt(N,-1,1);
tg = acos(xg);

fg = xg.^ltest - xg.^2;
sum_x_lmax = sum(wg.*fg)

% calcule les Pl0
plg = zeros(N,lmax+1);
pl = zeros(N,lmax+1);
for l=m:lmax
  ar=legendre(l,xg,'norm');
  plg(:,l+1) = ar(m+1,:)';
  ar=legendre(l,x,'norm');
  pl(:,l+1) = ar(m+1,:)';
end

% calcule zg
zg = plg .* (wg*ones(1,lmax+1));

% projette sur Plmax
for l=m:lmax
  cl(l+1) = sum(fg.*zg(:,l+1));
end

% synthetise
fgs = zeros(N,1);
fs = zeros(N,1);
for l=m:lmax
  fgs = fgs + cl(l+1) * plg(:,l+1);
  fs = fs + cl(l+1) * pl(:,l+1);
end
err_gauss = max(abs(fgs-fg))
err_dct = max(abs(fs-f))

figure; plot(tg,fg,'bo',tg,fgs,'g-');
figure; plot(t,f,'bo',t,fs,'g-');
figure; plot(t,fs-f);


%%% PROJECTION SUR GRILLE REGULIERE AVEC DCT + POINTS DE GAUSS %%%
for l=m:lmax
  for k=0:(N-1)
    z(k+1,l+1) = sum(zg(:,l+1) .* cos(k*tg))/N;
  end
  z(1,l+1) = z(1,l+1)/2;  % special case k=0
end

for l=m:lmax
    c2(l+1) = sum(df.*z(:,l+1));
end
errmax = max(abs(c2-cl))

%% Idem, avec idct_fftw ? %%%
z = z*4;
z(1,:) = z(1,:)/2;	% k=0 special case
zi = idct_fftw(z);
for l=m:lmax
    c3(l+1) = sum(f.*zi(:,l+1));
end
errmax = max(abs(c3-cl))


%%% BOUCLE avec donnees RANDOM %%%

f = -1+2*rand(size(t));		% random values
f0 = f;		% save random values

% FILTER
  % analyse
  for l=m:lmax
      cb(l+1) = sum(f.*zi(:,l+1));
  end
  % synthese
  f = zeros(N,1);
  for l=m:lmax
    f = f + cb(l+1)*pl(:,l+1);
  end
  f0 = f;

for ii=1:iter
  % analyse
  for l=m:lmax
      cb(l+1) = sum(f.*zi(:,l+1));
  end
  % synthese
  f = zeros(N,1);
  for l=m:lmax
    f = f + cb(l+1)*pl(:,l+1);
  end
  err(ii) = max(abs(f-f0));
end

figure; plot(t,f0,t,f);
figure; plot(log10(1:iter),log10(err));


%%% IDEM AVEC POINTS DE GAUSS %%%

f = -1+2*rand(size(t));		% random values
f0 = f;		% save random values

% FILTER
  % analyse
  for l=m:lmax
      cb(l+1) = sum(f.*zg(:,l+1));
  end
  % synthese
  f = zeros(N,1);
  for l=m:lmax
    f = f + cb(l+1)*plg(:,l+1);
  end
  f0 = f;

for ii=1:iter
  % analyse
  for l=m:lmax
      cb(l+1) = sum(f.*zg(:,l+1));
  end
  % synthese
  f = zeros(N,1);
  for l=m:lmax
    f = f + cb(l+1)*plg(:,l+1);
  end
  errg(ii) = max(abs(f-f0));
end

figure; plot(tg,f0,tg,f);
figure; plot(log10(1:iter),log10(err),log10(1:iter),log10(errg));
