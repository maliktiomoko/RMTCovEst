function [omega_mle,sigma_mle] = Lin_reg( X,Y,plot1 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
pt=1000;
x=linspace(min(X(:,1)),max(X(:,1)),pt);
X_norm=[X ones(size(X,1),1)];
omega_mle=(X_norm'*X_norm)^(-1)*X_norm'*Y;
sigma_mle=(1./(size(X,1)))*sum((Y-(omega_mle'*X_norm')').^2);
if plot1
figure
b=(-omega_mle(3)+0.5)
y_mle=-(omega_mle(1)./omega_mle(2))*x+(0.5-omega_mle(3))./omega_mle(2);
plot(x,y_mle,'r-')
hold on
plot(X(Y==0,1),X(Y==0,2),'bo')
plot(X(Y==1,1),X(Y==1,2),'g*')
xlabel('coordonée x1','Fontsize',15)
ylabel('coordonnée x2','Fontsize',15)
title('Classification des données du train A par la régression linéaire','Fontsize',15)
h=legend('proba=0.5','label 0',' label 1');
set(h,'Fontsize',15)
set(gca,'FOntsize',15)
x_li=max(X(:,1))-min(X(:,1));
y_li=max(X(:,2))-min(X(:,2));
axis([min(X(:,1))-x_li/10 max(X(:,1))+x_li/10 min(X(:,2))-y_li/10 max(X(:,2))+y_li/10])
hold off
%figure
%syms M N O
%f=(1/(sqrt(2*pi*sigma_mle)))*exp(-((1-omega_mle'*[M;N;1]).^2)/(2*sigma_mle))-0.5;
%ezplot(f);
%hold on
%plot(X(Y==0,1),X(Y==0,2),'bo')
%plot(X(Y==1,1),X(Y==1,2),'g*')
end
end

