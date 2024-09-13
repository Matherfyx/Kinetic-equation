function y = GetvgFlux(C0,D0,B0,qu,qv,h,N)
%% 周期边界
%compute g1+(g1R),g1-(g1L),g2+(g2R),g2-(g2L)
for i = 1:qv+1
   gi = Legendre(i-1,0,h);
   basisL(i) = gi(0);
   basisR(i) = gi(h);
end
g1R = sum(D0(:,1:qv+1).*basisL,2);
g1R(N+1) = g1R(1);
g1L = sum(D0(:,1:qv+1).*basisR,2);
g1L = [g1L(N);g1L];
g2R = sum(B0(:,1:qv+1).*basisL,2);
g2R(N+1) = g2R(1);
g2L = sum(B0(:,1:qv+1).*basisR,2);
g2L = [g2L(N);g2L];
%% Numerical Flux
% The vg- and vg+ flux,we note vgFulx(N+1,2)
vgFlux(:,1) = g1L;
vgFlux(:,2) = -g2R;
y = vgFlux;

% %%常数边界
% %compute g1+(g1R),g1-(g1L),g2+(g2R),g2-(g2L)
% for i = 1:qv+1
%    gi = Legendre(i-1,0,h);
%    basisL(i) = gi(0);
%    basisR(i) = gi(h);
% end
% g1R = sum(D0(:,1:qv+1).*basisL,2);
% g1R(N+1) = g1R(N);
% g1L = sum(D0(:,1:qv+1).*basisR,2);
% g1L = [g1L(1);g1L];
% g2R = sum(B0(:,1:qv+1).*basisL,2);
% g2R(N+1) = g2R(N);
% g2L = sum(B0(:,1:qv+1).*basisR,2);
% g2L = [g2L(1);g2L];
% %% Numerical Flux
% % The vg- and vg+ flux,we note vgFulx(N+1,2)
% vgFlux(:,1) = g1L;
% vgFlux(:,2) = -g2R;
% y = vgFlux;
end