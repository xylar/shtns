#!/usr/bin/octave

job='dj_ref_noj0xb0'
z=[0.4 0.7];


id = 0;

for i=1:150
  for j=1:length(z)
    fn = sprintf('./xspp poltorU_%04d.%s line 1000 0.001,0,%f 0.001,0,0',i,job,z(j));
    disp( fn )
    system(fn);
    a = load('o_line.0');
    s = a(:,1);
    wa(:,j) = a(:,8)./s;
  end

  da(2*id+1,:) = wa(:,1)';
  da(2*id+2,:) = wa(:,2)';
  id = id+1;
end;

save('-ascii',sprintf('mv_line.%s.txt',job),'da');

