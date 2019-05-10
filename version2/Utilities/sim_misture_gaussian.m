function [ simulated_var,col ] = sim_misture_gaussian(n,mu,sigma,M,p )
po=size(sigma,1);
X=zeros(n,po);
for i=1:n
    for k=1:M
        X(k,:)=mvnrnd(mu(k,:),sigma(:,:,k));
    end
    U=rand();
    Y=1:length(p);
    Y=p(1);
    for l=2:length(p)
       Y=[Y;p(l)+sum(p(1:l-1))];
    end
k=min(find((U-Y<0)));
if isempty(k)
    k=1;
end
simulated_var(i,:)=X(k,:);
col(i)=k;
end
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


end

