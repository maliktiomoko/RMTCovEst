clear all
clc
close all
s=1;
p=2;n2=4;
altitude=1;
c2=p/n2;
lambda=[1;2];
m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
zeta=zeros(p,1);
min_lambda=min(lambda);
 max_lambda=max(lambda);        
 maxV = max_lambda*1.5;
 step=1e-6;
for i=1:length(lambda)
    zeta_p=lambda(i);
    if i==1
        zeta_m=0;
    else
        zeta_m=lambda(i-1);
    end
    while abs(zeta_p-zeta_m)>1e-7*abs(lambda(p)-lambda(1))
                        zeta_=(zeta_p+zeta_m)/2;
                        if m(zeta_)<0
                            zeta_m=zeta_;
                        else
                            zeta_p=zeta_;
                        end
    end
zeta(i)=(zeta_p+zeta_m)/2;
end
 minV=zeta(1)/2;
test=@(z) (z-kappa(p+1));
for i=1:length(lambda)
    test=@(z) test(z).*(z-kappa(i))./(z-zeta(i));
end
F1=@(z) (1/s).*log(-s.*m(z)+1);
F2=@(z) (-m(z)).*log(-s.*m(z)+1);
F3=@(z) m(z);
contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]); 

kappa=zeros(p+1,1);lambda1=[0;lambda];
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if m(kappa_)>1/s
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa(1)=(kappa_p+kappa_m)/2;
for i=2:length(lambda1)
    kappa_p=lambda1(i);
    kappa_m=lambda1(i-1);
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if m(kappa_)>1/s
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa(i)=(kappa_p+kappa_m)/2;
end


real_integral1=(1/(2*pi*1j))*trapz(contour,F1(contour));
real_integral1_compare=sum(lambda-kappa(2:end));
real_integral2=(1/(2*pi*1j))*trapz(contour,F2(contour));
real_integral3=(1/(2*pi*1j))*trapz(contour,F3(contour));
real_integral=-(c2/p)*sum(log(lambda))+(c2/p)*sum(log(lambda-kappa(1)))-...
    (1-c2)*sum(log(kappa(2:end)./lambda))-c2+sum(lambda-kappa(2:end));