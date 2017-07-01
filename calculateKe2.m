clear
pkg load symbolic

syms z ln E I m c k Ke Me Ce x Qe P(t) Ke2

N=[1-3*z^2/ln^2+2*z^3/ln^3 z-2*z^2/ln+z^3/ln^2 3*z^2/ln^2-2*z^3/ln^3 -z^2/ln+z^3/ln^2];
Nx=[1-3*x^2/ln^2+2*x^3/ln^3 x-2*x^2/ln+x^3/ln^2 3*x^2/ln^2-2*x^3/ln^3 -x^2/ln+x^3/ln^2];
dN=diff(N,'z',1);
ddN=diff(N,'z',2);


% Ke2 : remove spring;
for i=1:4;
  for j=1:4;
    Ke2(i,j)=int(E*I*ddN(i)*ddN(j),'z',0,ln);
  end;
end;