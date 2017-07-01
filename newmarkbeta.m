% Newmark beta method.
% Reference: zeng qing yuan, liu jing bo;
% Powered by Duan Xiaoxu;
% K,M from case2_3 whose both ends' 2 vertical dofs is restrained;

clear;
syms t;
dofs=15;
deltaT=10;
delta=1/2;
alpha=1/4;
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

b0=1/(alpha*deltaT^2);
b1=delta/(alpha*deltaT);
b2=1/(alpha*deltaT);
b3=1/(2*alpha)-1;
b4=delta/alpha-1;
b5=deltaT/2*(delta/alpha-2);

equivalentK=K+b0*M;

while 1
   iteration=iteration+1;
   runtimeFlag='Caculating:'
   iteration
   t=t+deltaT;
   equivalentQt=double(eval(Qt))+M*(b0*qt+b2*Dqt+b3*D2qt);
   preqt=qt;
   preDqt=Dqt;
   preD2qt=D2qt;
   qt=inv(equivalentK)*equivalentQt;
   
   D2qt=b0*(qt-preqt)-b2*preDqt-b3*preD2qt;
   Dqt=b1*(qt-preqt)-b4*preDqt-b5*preD2qt;
   
   resultq=[resultq,qt];
   runtimeFlag='This iteration is completed.'
   runtimeFlag='------------------------------'
   
   %stop condition 1;
   if(iteration==stopT/deltaT)
      stopFlag='Stop at:'
      time=stopT;
      time
      break;    
   end    
end