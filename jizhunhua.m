% ji zhun hua matrix.
% used in subspace iteration.
function [res] = jizhunhua( matrix )
%JIZHUNHUA Summary of this function goes here
%  Detailed explanation goes here
    [row, col]=size(matrix);
    for i=1: col
        maxtemp=max(matrix(:,i));
        res(:,i)=matrix(:,i)/maxtemp;
    end    
end