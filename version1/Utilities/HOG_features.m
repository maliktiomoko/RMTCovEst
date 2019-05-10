function [features] = HOG_features(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
features=[];
for i=1:size(X,1)
    A=X(i,:);
    I=reshape(A,sqrt(784),sqrt(784));
    features = [features;extractHOGFeatures(I)];
end
end

