% Powered by Duan Xiaoxu;
% K,M from case2_3 whose ends are simple supported;
clear;

K=[
4374000000	0	-6561000	2187000000	0	0	0	0	0	0	0	0	0	0	0	;
0	3240000	0	0	-1620000	0	0	0	0	0	0	0	0	0	0	;
-6561000	0	26244	0	0	-13122	6561000	0	0	0	0	0	0	0	0	;
2187000000	0	0	8748000000	0	-6561000	2187000000	0	0	0	0	0	0	0	0	;
0	-1620000	0	0	3240000	0	0	-1620000	0	0	0	0	0	0	0	;
0	0	-13122	-6561000	0	26244	0	0	-13122	6561000	0	0	0	0	0	;
0	0	6561000	2187000000	0	0	8748000000	0	-6561000	2187000000	0	0	0	0	0	;
0	0	0	0	-1620000	0	0	3240000	0	0	-1620000	0	0	0	0	;
0	0	0	0	0	-13122	-6561000	0	26244	0	0	-13122	6561000	0	0	;
0	0	0	0	0	6561000	2187000000	0	0	8748000000	0	-6561000	2187000000	0	0	;
0	0	0	0	0	0	0	-1620000	0	0	3240000	0	0	-1620000	0	;
0	0	0	0	0	0	0	0	-13122	-6561000	0	26244	0	0	6561000	;
0	0	0	0	0	0	0	0	6561000	2187000000	0	0	8748000000	0	2187000000	;
0	0	0	0	0	0	0	0	0	0	-1620000	0	0	1620000	0	;
0	0	0	0	0	0	0	0	0	0	0	6561000	2187000000	0	4374000000	;

];

M=[
666666.6667	0	2166.6667	-500000	0	0	0	0	0	0	0	0	0	0	0	;
0	46.6667	0	0	11.6667	0	0	0	0	0	0	0	0	0	0	;
2166.6667	0	52	0	0	9	-2166.6667	0	0	0	0	0	0	0	0	;
-500000	0	0	1333333.333	0	2166.6667	-500000	0	0	0	0	0	0	0	0	;
0	11.6667	0	0	46.6667	0	0	11.6667	0	0	0	0	0	0	0	;
0	0	9	2166.6667	0	52	0	0	9	-2166.6667	0	0	0	0	0	;
0	0	-2166.6667	-500000	0	0	1333333.333	0	2166.6667	-500000	0	0	0	0	0	;
0	0	0	0	11.6667	0	0	46.6667	0	0	11.6667	0	0	0	0	;
0	0	0	0	0	9	2166.6667	0	52	0	0	9	-2166.6667	0	0	;
0	0	0	0	0	-2166.6667	-500000	0	0	1333333.333	0	2166.6667	-500000	0	0	;
0	0	0	0	0	0	0	11.6667	0	0	46.6667	0	0	11.6667	0	;
0	0	0	0	0	0	0	0	9	2166.6667	0	52	0	0	-2166.6667	;
0	0	0	0	0	0	0	0	-2166.6667	-500000	0	0	1333333.333	0	-500000	;
0	0	0	0	0	0	0	0	0	0	11.6667	0	0	23.3333	0	;
0	0	0	0	0	0	0	0	0	0	0	-2166.6667	-500000	0	666666.6667	;

];


% Subspace iteration to get the first r frequency;
nFreedoms=15;
r=5;
% frequency tolerance;
wTolerance=10^(-20);

% initialize psi0;
iteration=0;
psi=eye(nFreedoms, r)+1;
A=psi;
runtimeFlag='Initialized iteration 0;'
runtimeFlag='-----------'
% iteration;
while 1
    % record iterations;
    iteration=iteration+1;
    runtimeFlag='Calculating iteration: '
    iteration
    psi=inv(K)*M*A;
    psi=jizhunhua(psi);
    generalK=psi'*K*psi;
    generalM=psi'*M*psi;

    [a, lambda]=eig(inv(generalM)*generalK);
    % mode of vibration matrix;
    A=psi*a;
    A=jizhunhua(A);
    % frequency matrix;
    w=lambda.^(1/2);
    
    % record previous iteration wVector for tolerance check;
    if iteration~=1
        preWVector=wVector;   
    end    
    for i=1:r
        wVector(i,1)=w(i,i);
    end    
    % frequency vector from lower to higher;
    wVector=sort(wVector); 
    
    runtimeFlag='Iteration: '
    iteration
    runtimeFlag='completed;'
    runtimeFlag='-----------'
        
    % stop condition 1;
    if iteration ==1000
        stopFlag='Stop at iteration = 1000;'
        break;
    end
    % stop condition 2;
    if iteration~=1
        flag=0;
        for i=1:r
            if abs(wVector(i)-preWVector(i))/wVector(i)<=wTolerance
                flag=flag+1;    
            end    
        end 
        if flag==r
            stopFlag='Stop for sufficient w tolerance;'
            break;        
        end 
    end   
end 

% Compare with theoretical values;
errorW=checkError(wVector);
