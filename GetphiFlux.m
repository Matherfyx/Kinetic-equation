function y = GetphiFlux(C0,D0,B0,qu,qv,h,N)
%周期边界
%compute g1+(g1R),g1-(g1L),g2+(g2R),g2-(g2L)
for i = 1:qv+1
   phi = Legendre(i-1,0,h);
   basisL(i) = phi(0);
   basisR(i) = phi(h);
end
phiR= sum(C0(:,1:qv+1).*basisL,2);
phiR(N+1) = phiR(1);
phiL = sum(C0(:,1:qv+1).*basisR,2);
phiL = [phiL(N);phiL];
%% Numerical Flux
% The phi- flux,we note phiflux(N+1,1)
phiflux = phiL;
y = phiflux;

% %常数边界
% %compute g1+(g1R),g1-(g1L),g2+(g2R),g2-(g2L)
% for i = 1:qv+1
%    phi = Legendre(i-1,0,h);
%    basisL(i) = phi(0);
%    basisR(i) = phi(h);
% end
% phiR= sum(C0(:,1:qv+1).*basisL,2);
% phiR(N+1) = phiR(N);
% phiL = sum(C0(:,1:qv+1).*basisR,2);
% phiL = [phiL(1);phiL];
% %% Numerical Flux
% % The phi- flux,we note phiflux(N+1,1)
% phiflux = phiL;
% y = phiflux;
end