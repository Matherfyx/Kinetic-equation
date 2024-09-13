%% DG2-IMEX2时间离散
function [C,D,B] = DG_IMEX_RK2(C0,D0,B0,h,N,vgmflux,vgFlux,phiflux,dt,invA1,B1,basisMatrix,eps,qu,qv)
invA2 = invA1;
invA3 = invA1;
B2 = B1;
B3 = B1;
%% IMEXRK table
gamma = 1 - 1/sqrt(2);
delta = 1 - 1/(2*gamma);
ERK2 = [0,0,0;gamma,0,0;delta,1-delta,0];
Eb = [delta,1-delta,0];
Ec = [0;gamma;1];
IRK2 = [0,0,0;0,gamma,0;0,1-gamma,gamma];
Ib = [0,1-gamma,gamma];
Ic = [0;gamma;1];

%%
% SC(1,:,:) = C0;
% SD(1,:,:) = D0;
% SB(1,:,:) = B0;
% Svgmflux(1,:) = vgmflux;
% SvgFlux(1,:,:) = vgFlux;
% Sphiflux(1,:) = phiflux;
% Saux1 = zeros(3,N,qv+1);
% Saux2 = zeros(3,N,qv+1);
% csum_explicit = zeros(3,N,qu+1);

%% first stage
%% first we solve the first equation and solve the cofficient C for  the n+1 time
for i = 1:N
    SC(i,:) = C0(i,:)+ERK2(2,1)*dt*(invA1*(1./2*B1*D0(i,:)'-1./2*B1*B0(i,:)'-basisMatrix*[vgmflux(i+1);vgmflux(i)]))';
end
%% second we solve the auxiliary variable
for i = 1:N
    Saux1(i,:) = (invA2*(-B2*D0(i,:)'+basisMatrix*[vgFlux(i+1,1);vgFlux(i,1)]))';
    Saux2(i,:) = (invA2*(B2*B0(i,:)'+basisMatrix*[vgFlux(i+1,2);vgFlux(i,2)]))';
end
%% Third we solve the second equation and solve the cofficient D for the n+1 time
Sphiflux = GetphiFlux(SC,D0,B0,qu,qv,h,N);
for i = 1:N
    SD(i,:) = 1./(eps^2+IRK2(2,2)*dt)*(eps^2*D0(i,:)'-dt*(IRK2(2,1)*D0(i,:)'+eps/2*ERK2(2,1)*(Saux1(i,:)'-Saux2(i,:)')-invA3*B3*(IRK2(2,1)*C0(i,:)'+IRK2(2,2)*SC(i,:)')+invA3*basisMatrix*(IRK2(2,1)*[phiflux(i+1);phiflux(i)]+IRK2(2,2)*[Sphiflux(i+1);Sphiflux(i)])))';
    SB(i,:) = 1./(eps^2+IRK2(2,2)*dt)*(eps^2*B0(i,:)'-dt*(IRK2(2,1)*B0(i,:)'-eps/2*ERK2(2,1)*(Saux1(i,:)'-Saux2(i,:)')+invA3*B3*(IRK2(2,1)*C0(i,:)'+IRK2(2,2)*SC(i,:)')-invA3*basisMatrix*(IRK2(2,1)*[phiflux(i+1);phiflux(i)]+IRK2(2,2)*[Sphiflux(i+1);Sphiflux(i)])))';
end
Svgmflux = GetvgmFlux(SC,SD,SB,qu,qv,h,N);
SvgFlux = GetvgFlux(SC,SD,SB,qu,qv,h,N);


%% second stage
%% first we solve the first equation and solve the cofficient C for  the n+1 time
for i = 1:N
    a = ERK2(3,2)*dt*(invA1*(1./2*B1*SD(i,:)'-1./2*B1*SB(i,:)'-basisMatrix*[Svgmflux(i+1);Svgmflux(i)]))';
    C(i,:) = C0(i,:)+ a + ERK2(3,1)*dt*(invA1*(1./2*B1*D0(i,:)'-1./2*B1*B0(i,:)'-basisMatrix*[vgmflux(i+1);vgmflux(i)]))';
end
%% second we solve the auxiliary variable
for i = 1:N
    SSaux1(i,:) = (invA2*(-B2*SD(i,:)'+basisMatrix*[SvgFlux(i+1,1);SvgFlux(i,1)]))';
    SSaux2(i,:) = (invA2*(B2*SB(i,:)'+basisMatrix*[SvgFlux(i+1,2);SvgFlux(i,2)]))';
end
%C = C0 + dt*invA1*(reshape((1./2*B1*reshape(D0',[2*N,1])-1./2*B1*reshape(B0',[2*N,1])),[2,N])'-)
%C(:,1:qu+1) = C0(:,1:qu+1) + dt*(1./2*B1*D0(:,1:qv+1)'-1./2*B1*B0(:,1:qv+1)'-basisMatrix*fluxMatrix(:,1:(N+1)))';
%% Third we solve the second equation and solve the cofficient D for the n+1 time
SSphiflux = GetphiFlux(C,D0,B0,qu,qv,h,N);
for i = 1:N
    D(i,:) = 1./(eps^2+IRK2(3,3)*dt)*(eps^2*D0(i,:)'-dt*(IRK2(3,1)*D0(i,:)'+IRK2(3,2)*SD(i,:)'+eps/2*(ERK2(3,1)*(Saux1(i,:)'-Saux2(i,:)')+ERK2(3,2)*(SSaux1(i,:)'-SSaux2(i,:)'))-invA3*B3*(IRK2(3,1)*C0(i,:)'+IRK2(3,2)*SC(i,:)'+IRK2(3,3)*C(i,:)')+invA3*basisMatrix*(IRK2(3,1)*[phiflux(i+1);phiflux(i)]+IRK2(3,2)*[Sphiflux(i+1);Sphiflux(i)]+IRK2(3,3)*[SSphiflux(i+1);SSphiflux(i)])))';
    B(i,:) = 1./(eps^2+IRK2(3,3)*dt)*(eps^2*B0(i,:)'-dt*(IRK2(3,1)*B0(i,:)'+IRK2(3,2)*SB(i,:)'-eps/2*(ERK2(3,1)*(Saux1(i,:)'-Saux2(i,:)')+ERK2(3,2)*(SSaux1(i,:)'-SSaux2(i,:)'))+invA3*B3*(IRK2(3,1)*C0(i,:)'+IRK2(3,2)*SC(i,:)'+IRK2(3,3)*C(i,:)')-invA3*basisMatrix*(IRK2(3,1)*[phiflux(i+1);phiflux(i)]+IRK2(3,2)*[Sphiflux(i+1);Sphiflux(i)]+IRK2(3,3)*[SSphiflux(i+1);SSphiflux(i)])))';
end










% for s = 2:3
%     %% first we solve the first equation and solve the cofficient C for  the n+1 time
%     for i = 1:N
%         csum_explicit(s-1) = (invA1*(1./2*B1*SD(s-1,i,:)'-1./2*B1*SB(s-1,i,:)'-basisMatrix*[Svgmflux(s-1,i+1);Svgmflux(s-1,i)]))';
%         SC(s,i,:) = SC(s-1,i,:)+dt*(ERK2(s,1:s-1).*csum_explicit(1:s-1));
%     end
%     %% second we solve the auxiliary variable
%     for i = 1:N
%         Saux1(s-1,i,:) = (invA2*(-B2*SD(s-1,i,:)'+basisMatrix*[SvgFlux(s-1,i+1,1);SvgFlux(s-1,i,1)]))';
%         Saux2(s-1,i,:) = (invA2*(B2*SB(s-1,i,:)'+basisMatrix*[SvgFlux(s-1,i+1,2);SvgFlux(s-1,i,2)]))';
%     end
%     %% Third we solve the second equation and solve the cofficient D for the n+1 time
%     Sphiflux(s,:) = GetphiFlux(SC(s,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
%     for i = 1:N
%         SD(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*SD(s-1,i,:)'-dt*(sum(IRK2(s,1:s-1).*SD(1:s-1,i,:),2)+eps/2*(sum(ERK2(s,1:s-1).*(Saux1(1:s-1,i,:)-Saux2(1:s-1,i,:))',2))-invA3*B3*(sum(IRK2(s,1:s).*(SC(1:s,i,:)'),2))+invA3*basisMatrix*(sum(IRK2(s,1:s).*([Sphiflux(1:s,i+1);Sphiflux(1:s,i)]),2))))';
%         SD(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*SB(s-1,i,:)'-dt*(sum(IRK2(s,1:s-1).*SB(1:s-1,i,:),2)+eps/2*(-sum(ERK2(s,1:s-1).*(Saux1(1:s-1,i,:)-Saux2(1:s-1,i,:))',2))+invA3*B3*(sum(IRK2(s,1:s).*(SC(1:s,i,:)'),2))-invA3*basisMatrix*(sum(IRK2(s,1:s).*([Sphiflux(1:s,i+1);Sphiflux(1:s,i)]),2))))';
%         %SB(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*SB(s-1,i,:)'-dt*(-eps/2*(Saux1(s-1,i,:)'-Saux2(s-1,i,:)')+invA3*B3*SC(s,i,:)'-invA3*basisMatrix*[Sphiflux(s,i+1);Sphiflux(s,i)]))';
%     end
%     Svgmflux(s,:) = GetvgmFlux(SC(s-1,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
%     SvgFlux(s,:,:) = GetvgFlux(SC(s-1,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
% end


% for s = 2:3
%     %% first we solve the first equation and solve the cofficient C for  the n+1 time
%     for i = 1:N
%         csum_explicit(s-1,i,:) = (invA1*(1./2*B1*[SD(s-1,i,1);SD(s-1,i,2)]-1./2*B1*[SB(s-1,i,1);SB(s-1,i,2)]-basisMatrix*[Svgmflux(s-1,i+1);Svgmflux(s-1,i)]))';
%         SC(s,i,:) = SC(s-1,i,:)+dt*(ERK2(s,1:s-1).*csum_explicit(1:s-1,i,:));
%     end
%     %% second we solve the auxiliary variable
%     for i = 1:N
%         Saux1(s-1,i,:) = (invA2*(-B2*[SD(s-1,i,1);SD(s-1,i,2)]+basisMatrix*[SvgFlux(s-1,i+1,1);SvgFlux(s-1,i,1)]))';
%         Saux2(s-1,i,:) = (invA2*(B2*[SB(s-1,i,1);SB(s-1,i,2)]+basisMatrix*[SvgFlux(s-1,i+1,2);SvgFlux(s-1,i,2)]))';
%     end
%     %% Third we solve the second equation and solve the cofficient D for the n+1 time
%     Sphiflux(s,:) = GetphiFlux(SC(s,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
%     for i = 1:N
%         SD(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*[SD(s-1,i,1);SD(s-1,i,2)]-dt*(sum(IRK2(s,1:s-1).*SD(1:s-1,i,:),2)+eps/2*(sum(ERK2(s,1:s-1).*([Saux1(1:s-1,i,1)-Saux2(1:s-1,i,1);Saux1(1:s-1,i,2)-Saux2(1:s-1,i,2)]),2))-invA3*B3*(sum(IRK2(s,1:s).*([SC(1:s,i,1);SC(1:s,i,2)]),2))+invA3*basisMatrix*(sum(IRK2(s,1:s).*([Sphiflux(1:s,i+1);Sphiflux(1:s,i)]),2))))';
%         SD(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*[SB(s-1,i,1);SB(s-1,i,2)]-dt*(sum(IRK2(s,1:s-1).*SB(1:s-1,i,:),2)+eps/2*(-sum(ERK2(s,1:s-1).*([Saux1(1:s-1,i,1)-Saux2(1:s-1,i,1);Saux1(1:s-1,i,2)-Saux2(1:s-1,i,2)]),2))+invA3*B3*(sum(IRK2(s,1:s).*([SC(1:s,i,1);SC(1:s,i,2)]),2))-invA3*basisMatrix*(sum(IRK2(s,1:s).*([Sphiflux(1:s,i+1);Sphiflux(1:s,i)]),2))))';
%         %SB(s,i,:) = 1./(eps^2+IRK2(s,s)*dt)*(eps^2*SB(s-1,i,:)'-dt*(-eps/2*(Saux1(s-1,i,:)'-Saux2(s-1,i,:)')+invA3*B3*SC(s,i,:)'-invA3*basisMatrix*[Sphiflux(s,i+1);Sphiflux(s,i)]))';
%     end
%     Svgmflux(s,:) = GetvgmFlux(SC(s-1,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
%     SvgFlux(s,:,:) = GetvgFlux(SC(s-1,:,:),SD(s-1,:,:),SB(s-1,:,:),qu,qv,h,N);
% end
% C = SC(3,:,:);
% D = SD(3,:,:);
% B = SB(3,:,:);