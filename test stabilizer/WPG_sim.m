f=200;
time=[0:1/f:15]';

h=0.8;
g=9.81;
w=sqrt(g/h);

p_d=zeros(length(time),1);
p_d=[time p_d];

p_aux=zeros(length(time),1);
p_aux(5*f+1:10*f+1)=0.02;
p_aux=[time p_aux];

x_d=zeros(length(time),1);
x_d=[time x_d];

sx_d=zeros(length(time),1);
sx_d=[time sx_d];

q=[-13,-3,-w];
k1=-q(1)*q(2)-1;
k2=(q(1)+q(2))/w;
k3=w;