pho_L2_error = zeros(5,1);
pho_L2_order = zeros(4,1);
j_L2_error = zeros(5,1);
j_L2_order = zeros(4,1);
%计算pho的L2误差
for i = 1:5
    [pho_L2_error(i),j_L2_error(i),X,N,rhoNumerical,rho_exact] = KineticEqnDG2IMEX2(-pi,pi,1,2*(2^i),1,1,1e-6);
    i;
end
%计算pho的L2误差阶
for i = 2:5
    pho_L2_order(i-1) = log(pho_L2_error(i-1)/pho_L2_error(i))/log(2.0);
    j_L2_order(i-1) = log(j_L2_error(i-1)/j_L2_error(i))/log(2.0);
end
%绘图
plot(X(1:N),rhoNumerical(1:N),'oblue');
hold on
plot(X(1:N),rho_exact(1:N),'-red');