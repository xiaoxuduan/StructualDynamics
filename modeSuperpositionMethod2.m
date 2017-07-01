% Mode superposition method.
% Reference: zeng qing yuan, liu jing bo;
% Powered by Duan Xiaoxu;
% K,M from case2_3 whose both ends' 2 vertical dofs is restrained;

clear;

syms t tau;
resultq=[];

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

Qt=[
0	;
0	;
0	;
0	;
0	;
-6000*sin(0.02*pi*t)	;
0	;
0	;
0	;
0	;
0	;
0	;
0	;
0	;
0	;
    
];

% A is the mode of vibaration matrix;
% w is the frequency matrix;
B=inv(M)*K;
[A,w2]=eig(B);
w=w2.^(1/2);

% frequency vector;
vw=[];
for i=15:-1:1
    vw(15-i+1,1)=w(i,i);
end
vw=sort(vw);

E=200000;
I=5467500;
m=0.07;
L=5000;

% theoratical values from jie gou dong li xue jiang yi. zeng qing yuqn [M]. P66;
zengtw1=pi^2*sqrt(E*I/(m*L^4));
zengtw2=4*zengtw1;
zengtw3=9*zengtw1;
% error percent to jie gou dong li xue jiang yi. zeng qing yuqn [M]. P66;
zengerrorTw1=(vw(1)-zengtw1)/zengtw1;
zengerrorTw2=(vw(2)-zengtw2)/zengtw2;
zengerrorTw3=(vw(3)-zengtw3)/zengtw3;

% zheng jiao hua matrix;
% zeng qing yuan P55;
% zhengjiaoA: zheng jiao hua's A;
[row col]=size(w);
zhengjiaoA=A;
for i=1:col
   for j=i+1 : col
      if w(i,i)==w(j,j)
        c= -(zhengjiaoA(:,i)'*M*zhengjiaoA(:,i))/(zhengjiaoA(:,i)'*M*zhengjiaoA(:,j));    
        zhengjiaoA(:,j)=zhengjiaoA(:,i)+c*zhengjiaoA(:,j);
      end
   end    
end    

% create general M matrix, general K matrix;
[row col]=size(w);
for i=1:col
   generalM(i)=zhengjiaoA(:,i)'*M*zhengjiaoA(:,i);   
   generalK(i)=zhengjiaoA(:,i)'*K*zhengjiaoA(:,i);
end    

% create zheng ze mode of vibration;
% zhengjiaozeA: zheng jiao hua A, then zheng ze hua A here;
zhengjiaozeA=zhengjiaoA; 
[row col]=size(w);
for i=1:col
   zhengjiaozeA(:,i)=zhengjiaoA(:,i)/generalM(i);
end    

% Solve equation: D2
% Forced vibration without initial displacement and initial velocity;
% q=zhengjiaozeA*T;
[row col]=size(w);
Pt=zhengjiaozeA'*Qt;
for i=1:col
   if Pt(i,1)~=0
      % Duhamel integration;
      %T(i,1)=1/(generalM(i)*w(i,i))* (int(Pt(i,1)*sin(w(i,i)*(t-tau)),tau,0,t));
      T(i,1)=dsolve('D2TT+w(i,i)^2*TT-Pt(i,1)=0','TT(0)=0,DTT(0)=0','t');
   else
      T(i,1)=0;
   end    
end    
% superposion;
q=zhengjiaozeA*T;

% double result to resultq;
for t=0:10:300
   doubleq=double(eval(eval(q)));
   resultq=[resultq,doubleq];
end  

stopFlag='Program is over.'
