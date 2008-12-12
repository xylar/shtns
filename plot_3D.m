z0=-0.2;
fdisc = 'o_disc.0';	component=3;
fmerid = 'o_Vz.0';

figure;

% disc
a=load(fdisc);

s=size(a);
np=(s(2)-1)/3+1;
ip=0:(np-1);
ip(end)=0;
phi=2*pi*ip/(np-1);

r=a(:,1);
x=r*cos(phi);
y=r*sin(phi);
z=z0*ones(size(x));

vz=a(:,(component+1):3:end);
vz(:,end+1) = vz(:,1);

surf(x,y,z,-vz); shading interp;
hold;

% merid
b=load(fmerid);

ct = b(1,2:end);
st = sqrt(1-ct.*ct);
r = b(2:end,1);

x = r*st;
z = r*ct;
y = zeros(size(x));

vz = b(2:end,2:end);

surf(x,y,z,vz); shading interp;

% sphere wire
r=1;
theta = linspace(0,2*pi,600);
phi = 0;
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');	plot3(y,x,z,'k');

phi = pi/4;
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');	plot3(y,x,z,'k');

phi = pi/8;
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');	plot3(y,x,z,'k');

phi = 3*pi/8;
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');	plot3(y,x,z,'k');

theta = pi/2;
phi = linspace(0,2*pi,600);
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');

axis equal
theta = pi/4;
phi = linspace(0,2*pi,600);
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta)*ones(size(phi));
plot3(x,y,z,'k');	plot3(y,x,-z,'k');
