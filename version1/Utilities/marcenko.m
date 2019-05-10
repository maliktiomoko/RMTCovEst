function [ out ] = marcenko( z,c,lambda )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
m_av=1;
m=2;
while abs(m-m_av)>eps
    m_av=m;
    m=1./(lambda-c*lambda-z-c.*z.*m.*lambda);
end
out=m;
end

