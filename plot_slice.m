function plot_slice(S)

r = S(2:end,1);
ct = S(1,2:end);
st = sqrt(1.0-ct.*ct);

x = r*st;
y = r*ct;

pcolor(x,y,S(2:end,2:end));
%shading interp;
%axis equal off;

