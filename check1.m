clear qq xx ;
syms xx;
qq=zeros(15,31);
qq=qq+xx;
temp=0;
for t=0:10:300
   temp=temp+1;
   qq(:,temp)=eval(q) ;
end    
qq2=double(qq);

t=25;
q1=eval(q)