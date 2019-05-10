function [distance_mod,distance_real] = SCM_theo(p,n,x,C,f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Cs=x*x'/n;
metric=@(D) mean(f(eig(D\C)));
altitude=1;c=p/n;
    lambda=[(1-sqrt(c))^2;(1+sqrt(c))^2];
    min_lambda=min(lambda);
    max_lambda=max(lambda);        
    maxV = max_lambda*1.5;
    minV = min_lambda*.5;step=1e-4;
    contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]);        
    integrand=@(z) f(1./z).*marcenko(z,c,1);
    distance_mod=real((-1/(2*pi*1i))*(trapz(contour,integrand(contour))));
    distance_real=metric(Cs);
end

