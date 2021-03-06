% Wilson-theta method.
% Reference: zeng qing yuan, liu jing bo;
% Powered by Duan Xiaoxu;
% K,M from case2_3 whose both ends' 2 vertical dofs is restrained;

% something is wrong with it; May the book is wrong(zeng qing yuan's);

clear;
syms t;
dofs=15;
deltaT=10;
theta=1.4;
% stop at 300s;
stopT=300;
% resultq's ith column stores the q at t=0,deltaT,deltaT*2...
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
for i=dofs:-1:1
    vw(dofs-i+1,1)=w(i,i);
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

% iteration 0;
iteration=0;
t=0;
qt=zeros(dofs,1);
Dqt=zeros(dofs,1);
Q0=double(eval(Qt));
D2qt=inv(M)*(Q0-K*qt);
resultq=[resultq,qt];

a0=6/(theta*deltaT)^2;
a1=3/(theta*deltaT);
a2=2*a1;
a3=theta*deltaT/2;
a4=a0/theta;
a5=-a2/theta;
a6=1-3/theta;
a7=deltaT/2;
a8=deltaT^2/6;
equivalentK=K+a0*M;

while 1
   iteration=iteration+1;
   runtimeFlag='Caculating:'
   iteration
   t=t+theta*deltaT;
   %tTheta=theta*delta*iteration;
   equivalentQTheta=double(eval(Qt))+M*(a0*qt+a2*Dqt+2*D2qt);
   qtTheta=inv(equivalentK)*equivalentQTheta;
   
   preD2qt=D2qt;
   D2qt=a4*(qtTheta-qt)+a5*Dqt+a6*D2qt;
   Dqt=Dqt+a7*(D2qt+preD2qt);
   qt=qt+2*a7*Dqt+a8*(D2qt+2*preD2qt);
   
   resultq=[resultq,qt];
   runtimeFlag='This iteration is completed.'
   runtimeFlag='------------------------------'
   
   %stop condition 1;
   if(iteration==stopT/deltaT)
      stopFlag='Stop at:'
      stopT
      break;    
   end    
end