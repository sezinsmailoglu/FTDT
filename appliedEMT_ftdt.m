%Define the cell size ?z to use (in meters) 
dx = 0.015; %m
%Define the time step ?t to use (in seconds) 
dt = 2.5*10^(-11); %seconds
%Define how many cells will constitute the mesh
KE = 700;
KSource = 1;
%how many time steps for the simulation 
nsteps = 2000;
%dielectrics in the mesh
ex=zeros(1,KE);
hy=zeros(1,KE);
CA = zeros(1,KE);
CB = zeros(1,KE);
CC = zeros(1,KE);
sigma = zeros(1,KE);
er = ones(1,KE);
mur = ones(1,KE);
e0 = 8.854*10^(-12);
mu0 = 4*pi*10^(-7);

for k = 1 : KE
    if k >350
        er(k) = 8;
        sigma(k) = 0.15;
    end
end

for k =  1: KE
   CA(k) = ( 1 - sigma(k)*dt/(2*e0*er(k)))/(1+sigma(k)*dt/(2*e0*er(k)));
   CB(k) = ( 1 - dt/(dx*e0*er(k)))/(1+sigma(k)*dt/(2*e0*er(k)));
   CC(k) = (dt/(mu0*mur(k)*dx));
end

T = 10^(-9);
M=moviein(nsteps);
for t=0:dt:dt*nsteps
% E field loop
for k=2:KE-1
ex(k)=CA(k)*ex(k)+CB(k)*(hy(k-1)-hy(k));
end
% Source
if( t < 6*T)
ex(KSource)=exp(-5*((t-3*T)/T)^2);
else
    ex(KSource) = 0;
end
% H field loop
for k=1:KE-1
hy(k)=hy(k)+CC(k)*(ex(k)-ex(k+1));
end
plot(ex);

axis([1 KE -10^65 10^65]);
M(:) = getframe ;
end
movie(M,1);

