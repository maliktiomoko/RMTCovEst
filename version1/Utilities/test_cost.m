clear all
clc
close all
p_vec=[256;254;252;250;248;246;244;242;240];
n=260;
for i=1:length(p_vec)
    i
    p=p_vec(i);
c=p/n;
C=toeplitz(0.7.^(0:p-1));
M=toeplitz(0.1.^(0:p-1));
f=@(z) log(z).^2;
F=@(z) z.*(log(z).^2-2.*log(z)+2);
fg=zeros(p,n);
for k=1:n
     fg(:,k) = mvnrnd(zeros(1,p),C);
end
hatC  = fg*fg'/(n);
lambda=sort(eig(M\hatC));
m=@(z) c*mean(1./(lambda-z))-(1-c)./z;
mp=@(z) c*mean(1./((lambda-z).^2))+(1-c)./(z.^2);
lambda1=sort(eig(M*C));
m1=@(z) mean(1./(lambda1-z));
%integrand1=@(w) -f(w).*(m1(1./w)./(w.^2));
%integrand1=@(z) f(-m(z)).*mp(z).*(m1(-1./m(z))./(m(z).^2));
%integrand1=@(z) f(-m(z)).*(-z.*mp(z)/c+(1-1/c)*mp(z)./m(z));
integrand11=@(z) F(-m(z)).*(-1./c);
integrand12=@(z) (1/c)*f(-1./m(z)).*((1-c)*mp(z)./m(z));
integrand1=@(z) integrand11(z);
integrand2=@(z) integrand11(z)+integrand12(z);
integrand3=@(z) integrand12(z);
altitude=1;

            min_lambda=min(lambda);
            max_lambda=max(lambda);        
            maxV = max_lambda*1.5;
            minV = min_lambda*.5;step=1e-4;
            contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]);        
            integrand = @(z) F(-m(z)).*(m(z)+z.*mp(z));
            out=(real(trapz(contour,(1/(2*pi*1i*c))*integrand(contour))));
            out1(i)=(real(trapz(contour,(-1/(2*pi*1i))*integrand1(contour))));
            out2(i)=(real(trapz(contour,(-1/(2*pi*1i))*integrand2(contour))));
            out5(i)=(real(trapz(contour,(-1/(2*pi*1i))*integrand3(contour))));
 metric(i)=mean(log(eig(M\C)).^2);
end
 figure
 plot(abs(out1-metric),'r*-')
 hold on
 plot(abs(out2-metric),'go-')
 legend('old','new')
 figure
 plot(out5)