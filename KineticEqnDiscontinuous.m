%function [rho_L2_error,X,N,rhoNumerical,rho_exact] = KineticEqnDiscontinuous(xa,xb,t_end,N,qu,qv,eps)
%% DG1-IMEX1 Method to solve kinetic equation for summer school project
%% initial variable 
xa = -1;
xb = 1;
N = 160; 
h = (xb-xa)/N;
t = 0;
qu = 0;
qv = 0;
X = xa:h:xb;
X = X';
eps = 0.7;
dx = h;
dt = max(0.5*eps*dx,0.01*(dx^2));
t_end = 0.25;

%% initial condition
r = (-2)./(1+sqrt(1-4*(eps^2)));
phi = @(x) (-1<x & x<=0).*2.0 + ...
    (x>0 & x<1).*1.0;
g1 = @(x) 0.*x;  % v=1
g2 = @(x) 0.*x; % v=-1

%% initial coefficient matrix
C0 = zeros(N,qu+1);  % pho
D0 = zeros(N,qv+1);  % g1
B0 = zeros(N,qv+1);  % g2

%% initial projection
for p = 1:(qu+1)*N
    j = mod(p - 1,qu+1);
    i = (p - 1 - j)/(qu+1) + 1;
    fp = Legendre(j,X(i),X(i + 1)); % phi initial projection
    g = @(x) fp(x).*phi(x);
    F1(p) = integral(g,X(i),X(i + 1));
end
C0 = (reshape(F1,[qu+1,N]))';
for p = 1:(qu+1)*N
    j = mod(p-1,qv+1);
    i = (p-1-j)/(qv+1)+1;
    fp = Legendre(j,X(i),X(i+1));  % g1 g2 initial projection
    y1 = @(x) fp(x).*g1(x);
    y2 = @(x) fp(x).*g2(x);
    F2(p) = integral(y1,X(i),X(i+1));
    F3(p) = integral(y2,X(i),X(i+1));
end
D0 = (reshape(F2,[qv+1,N]))';
B0 = (reshape(F3,[qv+1,N]))';


vgmflux = GetvgmFlux(C0,D0,B0,qu,qv,h,N);
for i = 1:qu+1
   gi = Legendre(i-1,0,h);
   basisL(i) = gi(0);
   basisR(i) = gi(h);
end

%% the first eqn mass matrix and so on.
for j = 1:qu+1
    for i = 1:qu+1
        gi = Legendre(i-1,0,h);
        gj = Legendre(j-1,0,h);
        y = @(x) gi(x).*gj(x);
        A1(i,j) = integral(y,0,h);
    end
end
for i = 1:qu+1
    for j = 1:qv+1
        [gi,Dgi] = Legendre(i-1,0,h);
        [gj,Dgj] = Legendre(j-1,0,h);
        y = @(x) gj(x).*Dgi(x);
        B1(i,j) = integral(y,0,h);
    end
end
%% solve the auxiliary variable,get the  mass matrix and stiff matrix.
A2 = A1;
B2 = B1;

% as v = 1 and v = -1, get vgflux
vgFlux = GetvgFlux(C0,D0,B0,qu,qv,h,N);
basisMatrix = [basisR',-basisL'];
fluxMatrix = [vgmflux(2:N+1)',vgmflux(1);vgmflux'];
%E = speye(N);
invA1 = inv(A1);
invA2 = inv(A2);
aux1 = zeros(N,qv+1);
aux2 = zeros(N,qv+1);
%invA1 = kron(E,invA1);
%B1 = Kron(E,B1);
%% the second eqn mass matrix and so on.
A3 = A1;
B3 = B1;
invA3 = inv(A3);
phiflux = GetphiFlux(C0,D0,B0,qu,qv,h,N);


%% time loop
while(t<t_end)
    dt = min(dt,t_end-t);
    t = t + dt;
    %% first we solve the first equation and solve the cofficient C for  the n+1 time
    for i = 1:N
        C(i,:) = C0(i,:)+dt*(invA1*(1./2*B1*D0(i,:)'-1./2*B1*B0(i,:)'-basisMatrix*[vgmflux(i+1);vgmflux(i)]))';
    end
    %% second we solve the auxiliary variable 
    for i = 1:N
        aux1(i,:) = (invA2*(-B2*D0(i,:)'+basisMatrix*[vgFlux(i+1,1);vgFlux(i,1)]))';
        aux2(i,:) = (invA2*(B2*B0(i,:)'+basisMatrix*[vgFlux(i+1,2);vgFlux(i,2)]))';
    end
    %C = C0 + dt*invA1*(reshape((1./2*B1*reshape(D0',[2*N,1])-1./2*B1*reshape(B0',[2*N,1])),[2,N])'-)
    %C(:,1:qu+1) = C0(:,1:qu+1) + dt*(1./2*B1*D0(:,1:qv+1)'-1./2*B1*B0(:,1:qv+1)'-basisMatrix*fluxMatrix(:,1:(N+1)))';
    %% Third we solve the second equation and solve the cofficient D for the n+1 time
    phiflux = GetphiFlux(C,D0,B0,qu,qv,h,N);
    for i = 1:N
        D(i,:) = 1./(eps^2+dt)*(eps^2*D0(i,:)'-dt*(eps/2*(aux1(i,:)'-aux2(i,:)')-invA3*B3*C(i,:)'+invA3*basisMatrix*[phiflux(i+1);phiflux(i)]))';
        B(i,:) = 1./(eps^2+dt)*(eps^2*B0(i,:)'-dt*(-eps/2*(aux1(i,:)'-aux2(i,:)')+invA3*B3*C(i,:)'-invA3*basisMatrix*[phiflux(i+1);phiflux(i)]))';
    end
    C0 = C;
    D0 = D;
    B0 = B;
    vgmflux = GetvgmFlux(C0,D0,B0,qu,qv,h,N);
    vgFlux = GetvgFlux(C0,D0,B0,qu,qv,h,N);
end


%% compute rho or g1 left point value on element interval 
rhoNumerical = zeros(N,1);
for i=1:N
    sum0 = 0;
    for j=1:qu+1
        fi = Legendre(j-1,0,h);
        sum0 = sum0 + C0(i,j)*fi(0);
    end
    rhoNumerical(i) = sum0;
end
%% compute error
% pho_L2_error
rho_L2_error = 0;
for i = 1:N
    g = @(x) 0.*x;
    for j = 1:qu+1
        [fj,Dfj] = Legendre(j-1,X(i),X(i+1));
        g = @(x) g(x) + C0(i,j)*fj(x);
    end
    ff = @(x) 1./r*exp(r*t_end)*sin(x);
    h = @(x) (g(x) - ff(x)).*(g(x)-ff(x));
    rho_L2_error = rho_L2_error + integral(h,X(i),X(i+1));
end
rho_L2_error = sqrt(rho_L2_error);


%% plot
% rho exact solution
rho_exact = 1./r*exp(r*t_end).*sin(X);
% g1 exact solution
g1_exact = exp(r*t_end).*cos(X);
% g2 exact solution
g2_exact = -exp(r*t_end).*cos(X);
% plot(X,rho_L2_error);
 % plot(X(1:N),rho_exact(1:N),'red-');
 % hold on 
plot(X(1:N),rhoNumerical,'greeno')