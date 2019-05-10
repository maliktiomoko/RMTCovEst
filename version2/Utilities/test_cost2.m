clear all
clc
close all
p_vec=[240;220;200;180;160;140];
n=260;
f=@(t) log(t).^2;
%f=@(z) -(1/4)*log(z)+(1/2)*log(1+z)-(1/2)*log(2);
%f = @(t) (1/2)*t-(1/2)*log(t)-(1/2);
%f=@(t) t;
for i=1:length(p_vec)
    i
    p=p_vec(i);
c=p/n;
C=toeplitz(0.9.^(0:p-1));
M=toeplitz(0.1.^(0:p-1));
fg=zeros(p,n);
for k=1:n
     fg(:,k) = mvnrnd(zeros(1,p),C);
end
hatC  = fg*fg'/(n);
 metric(i)=mean(f(eig((M\C^(-1)))));
 out(i)=RMT_estim(M,hatC,n,'Inverse_Fisher');
end
 figure
 plot((metric).^2,'r*-')
 hold on
 plot((out),'b.-')