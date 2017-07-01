clear

pkg load symbolic
syms z ln E I m c k Ke Me Ce x Qe Pt Ke2 x l d A

N1=[1-x/l 0 0 x/l 0 0];
N2=[0 1-3*(x/l)^2+2*(x/l)^3 x*(1-x/l)^2 0 3*(x/l)^2-2*(x/l)^3 x*(x/l-1)*(x/l)];
N2d=[0 1-3*(d/l)^2+2*(d/l)^3 d*(1-d/l)^2 0 3*(d/l)^2-2*(d/l)^3 d*(d/l-1)*(d/l)];
N=[N1;N2];
% Nx=[1-3*x^2/ln^2+2*x^3/ln^3 x-2*x^2/ln+x^3/ln^2 3*x^2/ln^2-2*x^3/ln^3 -x^2/ln+x^3/ln^2];
dN2=diff(N2,'x',1);
ddN2=diff(N2,'x',2);
dN1=diff(N1,'x',1);

% Ke=int((E*I*ddN'*ddN+k*N'*N),'z',0,ln);
% Ke=zeros(4,4);
for i=1:6;
  for j=1:6;
    Ke(i,j)=int(E*(I*ddN2(i)*ddN2(j)+A*dN1(i)*dN1(j)),'x',0,l);
  end;
end;

% Me=int((m*N'*N),'z',0,ln);
% Me=zeros(4,4);
for i=1:6;
  for j=1:6;
    Me(i,j)=int((m*(N2(i)*N2(j)+N1(i)*N1(j))),'x',0,l);
  end;
end;

% Ce=int((c*N'*N),'z',0,ln);
% Me=zeros(4,4);
for i=1:6;
  for j=1:6;
    Ce(i,j)=int((c*(N2(i)*N2(j)+N1(i)*N1(j))),'x',0,l);
  end;
end;

% Qe
Qe=Pt*N2d';

