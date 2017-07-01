% K,M from case2_2 whose both ends' 2 vertical dofs is restrained;

clear;

K=[
4374000000	-6561000	2187000000	0	0	0	0	0	0	0	;
-6561000	26244	0	-13122	6561000	0	0	0	0	0	;
2187000000	0	8748000000	-6561000	2187000000	0	0	0	0	0	;
0	-13122	-6561000	26244	0	-13122	6561000	0	0	0	;
0	6561000	2187000000	0	8748000000	-6561000	2187000000	0	0	0	;
0	0	0	-13122	-6561000	26244	0	-13122	6561000	0	;
0	0	0	6561000	2187000000	0	8748000000	-6561000	2187000000	0	;
0	0	0	0	0	-13122	-6561000	26244	0	6561000	;
0	0	0	0	0	6561000	2187000000	0	8748000000	2187000000	;
0	0	0	0	0	0	0	6561000	2187000000	4374000000	;

];

M=[
666666.6667	2166.6667	-500000	0	0	0	0	0	0	0	;
2166.6667	52	0	9	-2166.6667	0	0	0	0	0	;
-500000	0	1333333.333	2166.6667	-500000	0	0	0	0	0	;
0	9	2166.6667	52	0	9	-2166.6667	0	0	0	;
0	-2166.6667	-500000	0	1333333.333	2166.6667	-500000	0	0	0	;
0	0	0	9	2166.6667	52	0	9	-2166.6667	0	;
0	0	0	-2166.6667	-500000	0	1333333.333	2166.6667	-500000	0	;
0	0	0	0	0	9	2166.6667	52	0	-2166.6667	;
0	0	0	0	0	-2166.6667	-500000	0	1333333.333	-500000	;
0	0	0	0	0	0	0	-2166.6667	-500000	666666.6667	;

];

Q=[
0	;
624	;
192000	;
5376	;
-768000	;
0	;
0	;
0	;
0	;
0	;

];

% v is the mode of vibaration matrix;
% w is the frequency matrix;
A=inv(M)*K;
[v,d]=eig(A);
w=d.^(1/2);

% frequency vector;
vw=[];
for i=8:-1:1;
    vw(8-i+1,1)=w(i,i);
end;
vw=sort(vw);

E=200000;
I=5467500;
m=0.07;
L=5000;

% theoratical values from ji suan jie gou li xue. zhu ci mian [M]2009. P170;
zhutw1=22.373/L^2*sqrt(E*I/m);
zhutw2=61.67/L^2*sqrt(E*I/m);

% error percent to ji suan jie gou li xue. zhu ci mian [M]2009. P170;
zhuerrorTw1=(vw(1)-zhutw1)/zhutw1;
zhuerrorTw2=(vw(2)-zhutw2)/zhutw2;

% theoratical values from jie gou dong li xue jiang yi. zeng qing yuqn [M]. P66;
zengtw1=pi^2*sqrt(E*I/(m*L^4));
zengtw2=4*zengtw1;
zengtw3=9*zengtw1;

% error percent to jie gou dong li xue jiang yi. zeng qing yuqn [M]. P66;
zengerrorTw1=(vw(1)-zengtw1)/zengtw1;
zengerrorTw2=(vw(2)-zengtw2)/zengtw2;
zengerrorTw3=(vw(3)-zengtw3)/zengtw3;



