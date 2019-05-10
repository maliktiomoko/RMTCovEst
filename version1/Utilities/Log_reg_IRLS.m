function [W] =Log_reg_IRLS( X,Y,plot1,lambda )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
pt=1000;
x=linspace(min(X(:,1)),max(X(:,1)),pt);

X_norm=[X ones(size(X,1),1)];
W=zeros(3,1);
eta=sigmoid(W'*X_norm');
eps=10^-8;
D=diag(eta.*(1-eta));
k=1;
W=[W W(:,k)+inv(X_norm'*D*X_norm)*X_norm'*(Y-eta')];
while abs(sum((W(:,k+1)-W(:,k))))>eps
    eta=sigmoid(W(:,k+1)'*X_norm');
    D=diag(eta.*(1-eta));
    W=[W W(:,k+1)+inv(X_norm'*D*X_norm-lambda)*X_norm'*((Y-eta')-lambda*sum(W(:,k+1)))];
    k=k+1;
end
if plot1
figure
y2=-(W(1,end)./W(2,end))*x-(W(3,end)/W(2,end));
plot(x,y2,'r-')
hold on
plot(X(Y==0,1),X(Y==0,2),'bo')
plot(X(Y==1,1),X(Y==1,2),'g*')
xlabel('coordonée x1','Fontsize',15)
ylabel('coordonnée x2','Fontsize',15)
title('Classification par la régression logistique du train A','Fontsize',15)
h=legend('proba=0.5','label 0',' label 1');
set(h,'Fontsize',15)
set(gca,'Fontsize',15)
x_li=max(X(:,1))-min(X(:,1));
y_li=max(X(:,2))-min(X(:,2));
axis([min(X(:,1))-x_li/10 max(X(:,1))+x_li/10 min(X(:,2))-y_li/10 max(X(:,2))+y_li/10])
end

end

