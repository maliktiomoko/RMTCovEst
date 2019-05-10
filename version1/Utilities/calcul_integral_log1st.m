clear all
clc
close all
 lambda=[2
    3];c2=0.5;p=length(lambda);
 m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;s=1;
 f=@(z) log(1+s*z);
 interior=@(z) f(-m(z));
%lambda=[3;4];zeta=lambda-0.5;
lambdak=lambda(1);
epsilon=1e-10;

integr=@(z) 1./((lambdak-z).^(2));
integral_real=0;
% Détermination des zeros de 1-s/m
kappa=zeros(p+1,1);lambdac=[0;lambda];
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa(1)=(kappa_p+kappa_m)/2;
for i=2:length(lambdac)
    kappa_p=lambdac(i);
    kappa_m=lambdac(i-1);
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa(i)=(kappa_p+kappa_m)/2;
end
%Détermination des zeros de m
zeta=zeros(p,1);
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
test=@(z) (z-kappa(p+1));
for i=1:length(lambda)
    test=@(z) test(z).*(z-kappa(i))./(z-zeta(i));
end
% figure
% z1=-2:0.01:10;
% plot(z1,interior(z1),'r*-');
% hold on
% plot(z1,test(z1),'go-')
% axis([min(z1) max(z1) -1000 1000])
lambdac=lambda;lambdac(1)=0;
integral_calcul=-sum(1./(zeta-lambdak)-1./(lambda(2:end)-lambdak));%+prod(lambdak-lambdac)./prod(lambdak-zeta)
altitude=1;
 min_lambda=min(lambda);
 max_lambda=max(lambda);        
 maxV = max_lambda*1.5;
 %minV = min_lambda*.2;
 minV=zeta(1)/2;
 step=1e-6;
contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]); 
integrand=@(z) (1/(2*pi*1j))*((interior(z)))./((lambdak-z).^(2));
real_integral=trapz(contour,integrand(contour))
calcul=prod(lambdak-zeta)
% compare1=0;compare2=0;compare3=0;eps=epsilon/(1e10);
% integral_arc1=0;integral_arc2=0;anc2=0;anc3=0;
% g11=@(theta) lambda(1)+epsilon*(cos(theta)+1j*sin(theta));
% der=@(theta) (g11(theta)-zeta(1))./g11(theta);
% for i=2:length(lambda)
%     der=@(theta) der(theta).*(g11(theta)-zeta(i))./(g11(theta)-lambda(i));
% end
% derm=@(x) (x-zeta(1))./x;
% for i=2:length(lambda)
%     derm=@(x) derm(x).*(x-zeta(i))./(x-lambda(i));
% end
% angle1=@(th) 1./(abs(lambda(1)-zeta(1)))-1./(abs(lambda(1)));
% for i=2:length(lambda)
%     angle1=@(th) angle1(theta)+1./(abs(lambda(1)-zeta(i)))-1./(abs(lambda(1)-lambda(i)));
% end
% % angle1=@(th) angle1(th).*epsilon.*sin(th);
% % figure
% % plot(angle(der(-pi+eps:pi/1000:pi-eps)))
% % hold on
% % plot(angle1((-pi+eps:pi/1000:pi-eps)))
% integral_arc_compare3=(1/(2*pi))*integral(@(theta)((log((1/epsilon).*der((theta)))).^2-theta.^2-2.*1j.*theta.*log((1/epsilon).*der((theta)))).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access1=(1/(2*pi))*integral(@(theta)((log((1/epsilon).*der((theta)))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access11=(1/(2*pi))*integral(@(theta)((log(abs((1/epsilon).*der((theta))))+1j*angle((1/epsilon)*der((theta)))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access2=(1/(2*pi))*integral(@(theta)((-theta.^2).*(1/epsilon).*exp(-1j*theta)),pi-eps,-pi+eps);
% first=@(th) (g11(th)-zeta(1))./g11(th);
% con=sum(1./(abs(lambda(1)-zeta)))-sum(1./abs((lambda(1)-lambda(2:end))))-1./(lambdak);
% co=(sum(1./(lambda(1)-zeta))-...
%         sum(1./(lambda(1)-lambda(2:end)))-...
%         1./(lambda(1)));
% %second=@(th) zeta(1)*epsilon*sin(th)./((lambda(1)-zeta(1))*lambda(1));
% access111=(1/(2*pi))*integral(@(theta)((log(abs((1/epsilon).*der((theta))))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access1111=(1/(2*pi))*integral(@(theta)((log(abs((1/epsilon).*der((0))))+1j*epsilon*co*sin(theta)).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% %access1112=(1/(2*pi))*integral(@(theta)((log(abs((1/epsilon).*der((0))))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% %access1113=(1/(2*pi))*integral(@(theta)((log(abs((1/epsilon).*der((0))))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access112=(1/(2*pi))*integral(@(theta)((1j*angle(der((theta)))).^2).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access113=(1/(2*pi))*integral(@(theta)((2.*log(abs((1/epsilon).*der((theta)))).*1j.*angle(der((theta))))).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access1131=(-2/(pi)).*log(abs((1/epsilon).*derm((lambda(1))))).*co*((pi-eps)-0.5*sin(2*eps));
%  figure
%  plot(co.*epsilon.*sin(-pi+eps:pi/1000:pi-eps));
%  hold on
%  plot(angle(der((-pi+eps:pi/1000:pi-eps))))
% access3=(1/(2*pi))*integral(@(theta)(-2.*1j.*theta.*log((1/epsilon).*der((theta)))).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% access31=(1/(2*pi))*integral(@(theta)(-2.*1j.*theta.*(log(abs(1/epsilon).*derm((lambda(1)))))).*(1/epsilon).*exp(-1j*theta),pi-eps,-pi+eps);
% acess=access1+access2+access3;
% integral_arc_compare2=(-2/(pi)).*log(abs((1/epsilon).*derm((lambda(1))))).*co*((pi-eps)-0.5*sin(2*eps))-...
%     (2/epsilon)+2*log(abs((1/epsilon)*derm(lambda(1))))*(1/epsilon);
% %compare2=compare2+prim1(lambdak-epsilon)-prim1(zeta(1)-epsilon);
% %2*integral_real
% %compare1+compare3
% %2*(compare2+compare3)
% integral_arc2;
% %2*(compare2+compare3)-integral_arc2;
% %compare3
% % A=0;integral_real_compare=0;
% % for i=1:length(lambda)
% %     if lambda(i)~=lambdak
% %         A=A+1/(lambda(i)-lambdak)+1/(lambdak-zeta(i));
% %     end
% % end
% % A=A-1/lambdak;
% % for j=2
% %     if lambda(j)~=lambdak
% %     integral_real_compare=integral_real_compare+A*log(abs((lambda(j)-epsilon-lambdak)/(zeta(j)+epsilon-lambdak)))-(((lambda(j)-epsilon)*log(abs(lambda(j)-epsilon)))/(lambdak*(lambdak-lambda(j)+epsilon))+((zeta(j)+epsilon)*log((zeta(j)+epsilon)))/(lambdak*(lambdak-zeta(j)-epsilon)));
% %     %integral_real_compare=integral_real_compare+((zeta(1)-lambdak)*(log(abs((lambda(j)-epsilon-zeta(1))/(lambda(j)-lambdak-epsilon)))-1))/((zeta(1)-lambdak)*(lambdak-lambda(j)+epsilon));
% %     %integral_real_compare=integral_real_compare+((zeta(1)-zeta(j)-epsilon)*(log(abs((zeta(j)+epsilon-zeta(1))/(zeta(j)-lambdak-epsilon)))-1))/((zeta(1)-lambdak)*(lambdak-zeta(j)-epsilon));
% %     end
% % end
% %integral_real_compare
% % (1/(2*pi)).*(1/epsilon).*integral(@(theta)((log((1/epsilon).*der(g11(theta)))).^2).*exp(-1j*theta),pi-eps,-pi+eps);
% % (1/(2*pi)).*(1/epsilon).*integral(@(theta)(((log(1/epsilon).^2+log(der(g11(theta))).^2+2*log(1/epsilon).*log(der(g11(theta)))))).*exp(-1j*theta),pi-eps,-pi+eps);
% % int1=(1/(2*pi)).*(1/epsilon).*integral(@(theta)((log((1/epsilon)).^2).*exp(-1j*theta)),pi-eps,-pi+eps);
% % int2=(1/(2*pi)).*(1/epsilon).*integral(@(theta)((log(der(g11(lambda(1)))).^2).*exp(-1j*theta)),pi-eps,-pi+eps);
% % int3=(1/(pi)).*(1/epsilon).*integral(@(theta)((log((1/epsilon)).*log(der(g11(lambda(1))))).*exp(-1j*theta)),pi-eps,-pi+eps);
% % (1/(2*pi)).*((log((1/epsilon).*der(g11(lambda(1))))).^2).*(1/epsilon)*integral(@(theta) exp(-1j*theta),pi-eps,-pi+eps);
% real_integral
% %theta=-pi+eps:pi/1000:pi-eps;
% %theta=[-pi+eps pi/2];
% % x=lambda(1)+epsilon*exp(1j*theta);
% % der2=@(x) prod(x-zeta)./(x.*prod(x-lambda(2:end)));
% % der3=@(x) prod(x-zeta)./(x.*(x-lambda(2:end)));
% % plot(theta,imag(der2(x)),'go')
% % hold on
% % plot(theta,imag(der3(x)),'r*')
% % hold on
% % plot(theta,imag(der(x)),'k.')
% % second=@(th) (epsilon*sum(1./(lambda(1)-zeta))*sin(th));
% % second=@(th) zeta(1)*epsilon*sin(th)./((lambda(1)-zeta(1))*lambda(1));
% % figure
% % plot(-pi+eps:pi/1000:pi-eps,angle(first(-pi+eps:pi/1000:pi+eps)),'r-')
% % hold on
% % plot(-pi+eps:pi/1000:pi+eps,angle(der(-pi+eps:pi/1000:pi+eps)),'g-')
% % plot(-pi+eps:pi/1000:pi+eps,second(-pi+eps:pi/1000:pi-eps),'c.')
% % hold on
% % plot(-pi+eps:pi/1000:pi+eps,second(-pi+eps:pi/1000:pi+eps),'go')
% mp=@(z) 1./(z.*(z-lambda(1)));
% for i=2:length(lambda)
%     mp=@(z) mp(z).*(z-zeta(i))./(z-lambda(i));
% end
% addi=0;
% for i=2:length(lambda)
%     for j=2:length(lambda)
%         addi=addi+log(abs(((lambda(j)-epsilon-zeta(i))/(lambda(j)-epsilon-lambda(i)))))/(lambdak-lambda(j)+epsilon)-...
%             log(abs(((zeta(j)+epsilon-zeta(i))/(zeta(j)+epsilon-lambda(i)))))/(lambdak-zeta(j)-epsilon)+...
%             log(abs(((lambda(j)-epsilon-zeta(i))/(zeta(j)+epsilon-zeta(i)))))/(zeta(i)-lambdak)+...
%             log(abs(((lambda(j)-epsilon-lambda(i))/(zeta(j)+epsilon-lambda(i)))))/(lambdak-lambda(i));
%     end
% end
% addi=addi+sum(log(abs(((lambda(1)-epsilon-zeta(2:end))./(lambda(1)-epsilon-lambda(2:end)))))./(epsilon)-...
%             log(abs(((zeta(1)-zeta(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-zeta(1))+...
%             log(abs(((lambda(1)-zeta(2:end))./(zeta(1)-zeta(2:end)))))./(zeta(2:end)-lambdak)+...
%             log(abs(((lambda(1)-lambda(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-lambda(2:end)));
% %conn2=sum((zeta(2:end)-lambda(2:end))./((lambda(1)-lambda(2:end)).*(lambda(1)-zeta(2:end))));
% verif2=sum(((zeta+epsilon).*log(zeta+epsilon))./(lambda(1)*(lambda(1)-zeta-epsilon)))-...
%     sum(((lambda(1:end)-epsilon).*log(lambda(1:end)-epsilon))./(lambda(1).*(lambda(1)-lambda(1:end)+epsilon)))-...
%     (1/lambda(1))*sum((log(abs((lambda(1:end)-lambda(1)-epsilon)./(zeta(1:end)+epsilon-lambda(1))))));
% conn=sum((zeta(2:end)-lambda(2:end))./((lambda(1)-lambda(2:end)).*(lambda(1)-zeta(2:end))));
% conn2=-1./lambda(1)+sum((zeta(2:end)-lambda(2:end))./((lambda(1)-lambda(2:end)).*(lambda(1)-zeta(2:end))));%+...
%     verif3=conn*sum(log(abs((lambda-epsilon-lambda(1))./(zeta+epsilon-lambda(1)))));
% % tot=2*(-sum(log(abs(((zeta(1)-zeta(2:end)-epsilon)./(zeta(1)-lambda(2:end)+epsilon))))./(lambdak-zeta(1)))+...
% %             sum(log(abs(((lambda(1)-zeta(2:end)-epsilon)./(zeta(1)-zeta(2:end)-epsilon))))./(zeta(2:end)-lambdak+epsilon))+...
% %             sum(log(abs(((lambda(1)-lambda(2:end)+epsilon)./(zeta(1)-lambda(2:end)+epsilon))))./(lambdak-lambda(2:end)+epsilon))+...
% %             sum(((zeta+epsilon).*log(zeta+epsilon))./(lambda(1)*(lambda(1)-zeta-epsilon)))-...
% %     sum(((lambda(2:end)-epsilon).*log(lambda(2:end)-epsilon))./(lambda(1).*(lambda(1)-lambda(2:end)+epsilon)))+...
% %     conn2*sum(log(abs((lambda(2:end)-epsilon-lambda(1))./(zeta(2:end)+epsilon-lambda(1)))))+...
% %     sum((((zeta(1)-lambda(2:end)+epsilon).*(log(abs((lambda(2:end)-epsilon-zeta(1))./(lambda(2:end)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(2:end)+epsilon))))-...
% %     sum((((zeta(1)-zeta(2:end)-epsilon).*(log(abs((zeta(2:end)+epsilon-zeta(1))./(zeta(2:end)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(2:end)-epsilon))))+...
% %     conn2*sum(log(abs((epsilon)./(zeta(1)+epsilon-lambda(1)))))+...
% %     sum(log(abs(((lambda(1)-epsilon-zeta(2:end))./(lambda(1)-epsilon-lambda(2:end)))))./(epsilon))-...
% %     log(lambda(1)-epsilon)./epsilon-...
% %     sum((((-epsilon).*(log(abs((+epsilon)./(zeta(1)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(1)-epsilon))))+...
% %     sum((((zeta(1)-lambda(1)+epsilon).*(log(abs((lambda(1)-epsilon-zeta(1))./(-epsilon)))-1))./((zeta(1)-lambda(1)).*(epsilon)))))-...
% %     ((-2/(pi)).*log(abs((1/epsilon).*der((0)))).*con*((pi-eps)-0.5*sin(2*eps))-...
% %     (2/epsilon)+2*log(abs((1/epsilon)*der(lambda(1))))*(1/epsilon));
% 
% 
% 
% % verif=sum(log(abs(((lambda(1)-epsilon-zeta(2:end))./(lambda(1)-epsilon-lambda(2:end)))))./(epsilon)-...
% %             log(abs(((zeta(1)-zeta(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-zeta(1))+...
% %             log(abs(((lambda(1)-zeta(2:end))./(zeta(1)-zeta(2:end)))))./(zeta(2:end)-lambdak)+...
% %             log(abs(((lambda(1)-lambda(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-lambda(2:end)))+...
% %             conn*sum(log(abs(((lambda-epsilon-lambdak)./(zeta+epsilon-lambdak)))));
% 
% % int_cal=2*(-sum(log((zeta(1)-zeta(2:end))./(zeta(1)-lambda(2:end))))/(lambda(1)-zeta(1))+...
% %     conn2*sum((log(lambda(2:end)-lambda(1))./(zeta(2:end)-lambda(1))))+...
% %     sum((zeta.*log(zeta))./(lambda(1)*(lambda(1)-zeta)))-...
% %     sum((lambda(2:end).*log(lambda(2:end)))./(lambda(1).*(lambda(1)-lambda(2:end))))+...
% %     sum(log((lambda(1)-zeta(2:end))./(zeta(1)-zeta(2:end)))./(zeta(2:end)-lambda(1)))+...
% %     sum(log((lambda(1)-lambda(2:end))./(zeta(1)-lambda(2:end)))./(-lambda(2:end)+lambda(1)))-...
% %     sum((zeta(1)-zeta(2:end)).*(log((zeta(2:end)-zeta(1))./(zeta(2:end)-lambda(1)))-1)./((zeta(1)-lambda(1)).*(lambda(1)-zeta(2:end))))+...
% %     sum((zeta(1)-lambda(2:end)).*(log((lambda(2:end)-zeta(1))./(lambda(2:end)-lambda(1)))-1)./((zeta(1)-lambda(1)).*(lambda(1)-lambda(2:end)))))+...
% %     sum(log((lambda(1)-zeta(2:end))./(lambda(1)-lambda(2:end))))/epsilon-...
% %     (log(lambda(1)))/epsilon+log(epsilon./(zeta(1)-lambda(1)))/lambda(1)+...
% %     (zeta(1)-lambda(1))*(log((lambda(1)-zeta(1))/(epsilon)))/((zeta(1)-lambda(1))*(lambda(1)-zeta(1)))
% %     sum(log((zeta(2:end)-lambda(1))./(zeta(2:end)-zeta(1))))/(zeta(1)-lambda(1))+...
% %     sum(log((lambda(2:end)-zeta(1))./(lambda(2:end)-lambda(1))))/(zeta(1)-lambda(1))-...
% %     sum(log((lambda(2:end)-lambda(1))./(zeta(2:end)-zeta(1))))/lambda(1)
% cer=sum((((zeta(1)-lambda(2:end)+epsilon).*(log(abs((lambda(2:end)-epsilon-zeta(1))./(lambda(2:end)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(2:end)+epsilon))))+...
%     sum((((zeta(1)-lambda(1)+epsilon).*(log(abs((lambda(1)-epsilon-zeta(1))./(lambda(1)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(1)+epsilon))))-...
%         sum((((zeta(1)-zeta(2:end)-epsilon).*(log(abs((zeta(2:end)+epsilon-zeta(1))./(zeta(2:end)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(2:end)-epsilon))))-...
%         sum((((zeta(1)-zeta(1)-epsilon).*(log(abs((zeta(1)+epsilon-zeta(1))./(zeta(1)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(1)-epsilon))));
% tot=2*(-sum(log(abs(((zeta(1)-zeta(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-zeta(1)))+...
%             sum(log(abs(((lambda(1)-zeta(2:end))./(zeta(1)-zeta(2:end)))))./(zeta(2:end)-lambdak))+...
%             sum(log(abs(((lambda(1)-lambda(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-lambda(2:end)))+...
%             0+...
%             sum(((zeta+epsilon).*log(zeta+epsilon))./(lambda(1)*(lambda(1)-zeta-epsilon)))-...
%     sum(((lambda(2:end)-epsilon).*log(lambda(2:end)-epsilon))./(lambda(1).*(lambda(1)-lambda(2:end)+epsilon)))-...
%     (1/lambda(1))*sum((log(abs((lambda(1:end)-lambda(1)-epsilon)./(zeta(1:end)+epsilon-lambda(1))))))+...
%     0+...
%     conn*sum(log(abs((lambda-epsilon-lambda(1))./(zeta+epsilon-lambda(1)))))+...
%     0+...
%     sum((((zeta(1)-lambda(2:end)+epsilon).*(log(abs((lambda(2:end)-epsilon-zeta(1))./(lambda(2:end)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(2:end)+epsilon))))-...
%         sum((((zeta(1)-zeta(2:end)-epsilon).*(log(abs((zeta(2:end)+epsilon-zeta(1))./(zeta(2:end)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(2:end)-epsilon)))))+...
%         0+...
%         ((2).*log(abs((1/epsilon).*derm((lambda(1))))).*co-2*log(abs((1/epsilon)*derm(lambda(1))))*(1/epsilon))+...
%     2*sum(log(abs(((lambda(1)-epsilon-zeta(2:end))./(lambda(1)-epsilon-lambda(2:end)))))./(epsilon))-...
%     2*sum(((lambda(1)-epsilon).*log(lambda(1)-epsilon))./(lambda(1).*(epsilon)))+...
%     2*((((1).*(log(abs((lambda(1)-epsilon-zeta(1))./(-epsilon)))))./(1.*(epsilon))))-2/(zeta(1)-lambda(1))+...
%     2*((((1/(zeta(1)-lambda(1))).*(log(abs((lambda(1)-epsilon-zeta(1))./(epsilon)))))))
%     %2*((((zeta(1)-lambda(1)+epsilon).*(log(abs((lambda(1)-epsilon-zeta(1))./(-epsilon)))))./((zeta(1)-lambda(1)).*(epsilon))))-2/(zeta(1)-lambda(1))
% %2*((((zeta(1)-lambda(1)+epsilon).*1)./((zeta(1)-lambda(1)).*(epsilon))))
%     %2*((((-epsilon).*(log(abs((epsilon)./(zeta(1)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(1)-epsilon))));
%     tot-real_integral
% integral_real-integral_arc2
% a0=2*(-sum(log(abs(((zeta(1)-zeta(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-zeta(1)))+...
%             sum(log(abs(((lambda(1)-zeta(2:end))./(zeta(1)-zeta(2:end)))))./(zeta(2:end)-lambdak))+...
%             sum(log(abs(((lambda(1)-lambda(2:end))./(zeta(1)-lambda(2:end)))))./(lambdak-lambda(2:end)))+...
%             0+...
%             sum(((zeta+epsilon).*log(zeta+epsilon))./(lambda(1)*(lambda(1)-zeta-epsilon)))-...
%     sum(((lambda(2:end)-epsilon).*log(lambda(2:end)-epsilon))./(lambda(1).*(lambda(1)-lambda(2:end)+epsilon)))-...
%     (1/lambda(1))*sum((log(abs((lambda(1:end)-lambda(1)-epsilon)./(zeta(1:end)+epsilon-lambda(1))))))+...
%     0+...
%     conn*sum(log(abs((lambda-epsilon-lambda(1))./(zeta+epsilon-lambda(1)))))+...
%     0+...
%     sum((((zeta(1)-lambda(2:end)+epsilon).*(log(abs((lambda(2:end)-epsilon-zeta(1))./(lambda(2:end)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(2:end)+epsilon))))-...
%         sum((((zeta(1)-zeta(2:end)-epsilon).*(log(abs((zeta(2:end)+epsilon-zeta(1))./(zeta(2:end)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(2:end)-epsilon)))))+...
%         0;
% a1=((-2/(pi)).*log(abs((1/epsilon).*derm((lambda(1))))).*co*((pi-eps)-0.5*sin(2*eps))-...
%     (2/epsilon)+2*log(abs((1/epsilon)*derm(lambda(1))))*(1/epsilon));
% a2=2*sum(log(abs(((lambda(1)-epsilon-zeta(2:end))./(lambda(1)-epsilon-lambda(2:end)))))./(epsilon))-...
%     2*sum(((lambda(1)-epsilon).*log(lambda(1)-epsilon))./(lambda(1).*(lambda(1)-lambda(1)+epsilon)))+...
%     2*sum((((zeta(1)-lambda(1)+epsilon).*(log(abs((lambda(1)-epsilon-zeta(1))./(lambda(1)-epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-lambda(1)+epsilon))))-...
%     2*sum((((zeta(1)-zeta(1)-epsilon).*(log(abs((zeta(1)+epsilon-zeta(1))./(zeta(1)+epsilon-lambda(1))))-1))./((zeta(1)-lambda(1)).*(lambdak-zeta(1)-epsilon))));
% p=2;n=300;c2=p/n;
% out_prim=zeros(p,p);
% for i=1:length(lambda)
%     lambdac=lambda;
%     lambdac(i)=[];
% conn3=-1./lambda(i)+sum(1./(lambda(i)-zeta(1:end)))-sum(1./(lambda(i)-lambdac));
% c=-log(1-c2);
% out_prim(i,i)=(-2*conn3)*(log(lambda(i)))+(2+2*c)/lambda(i)+...
%             2*sum((1./(lambda(i)-zeta)).*log(zeta))-2*sum((1./(lambda(i)-lambdac)).*log(lambdac));
% end
% der1=@(theta) (g11(theta)-zeta(1))./g11(theta);
% for i=2:length(lambda)
%     der1=@(theta) der1(theta).*(g11(theta)-zeta(i))./(g11(theta)-lambda(i));
% end
% the=-pi+eps:pi/1000:pi-eps;
% vec_test=angle((der1(the)));
% for i=1:length(the)
%     %vec_test2(i)=sum(angle(g11(the(i))-zeta));
%         vec_test2(i)=epsilon*sin(the(i))*co;
% end
% % figure
% % plot(angle(der(the)))
% % hold on
% % plot(vec_test2)
te=zeta(1)*(prod(zeta(1)-lambda)./prod(zeta(1)-zeta(2:end)))/(lambdak-zeta(1))^2+zeta(2)*prod(zeta(2)-lambda)./prod(zeta(2)-zeta(1))/((lambdak-zeta(2))^2)+...
    lambdak*(lambdak-lambda(2))./((lambdak-zeta(1)).*(lambdak-zeta(2)))