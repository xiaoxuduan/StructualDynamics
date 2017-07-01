% Comparing with theoretical values;
% used in subspace.
function [res] = checkError(wVector)
%CHECKERROR Summary of this function goes here
%  Detailed explanation goes here
    E=200000;
    I= 5467500;
    m=0.07;
    L=5000;
    [row col]=size(wVector);
    for i=1:row
        theoryW(i)=i^2*pi^2*(E*I/(m*L^4))^0.5;
        res(i)=(wVector(i,1)-theoryW(i))/theoryW(i);
    end    
end
