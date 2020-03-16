clc;
clear all;
syms x
load('WH_delta1=3_delta3=3_k=0.3_alpha=0.05.mat')
w=Win;
delta1=3;
delta3=3;
kk=0.3;
alpha=0.05;
tau1=0.5;
u=0.01;
v=0.01;
hx=81;
hy=81;
hd=0.049;
x=-2+(0:hx-1)*hd;
n=length(x);
Ntau=[0:0.1:6];
n1=length(Ntau);
U=zeros(n1,n); % the potential function
deltaU=zeros(n1,1); % the depth of potential wells
for j=1:n1
    tau2=Ntau(j);
    for i=1:n
        w1=w(i,42);
        U(j,i)=-1/2*delta1*x(i)^2+1/4*delta3*x(i)^4+1/2*(kk*w1^2/(alpha^2+w1^2)-u*cos(w1*tau1)-v*w1*sin(w1*tau2))*x(i)^2;
    end
    deltaU(j)=abs(min(U(j,:)));
end

figure;
plot(Ntau,deltaU,'r-*');%the vary of well depth
xlabel('$\tau_2$','Interpreter','latex');
ylabel('$\Delta U$','Interpreter','latex');

% figure;
% plot(x,U,'r-*'); %the vary of potential function
% xlabel('$X$','Interpreter','latex');
% ylabel('$U(X)$','Interpreter','latex');
% ylim([-0.8,0.1]);

figure;
mesh(Ntau,x(15:69),U(:,15:69)') %the vary of potential function
xlabel('$\tau_2$','Interpreter','latex');
ylabel('$X$','Interpreter','latex');
zlabel('$U(X)$','Interpreter','latex');
ylim([-1.5,1.5]);
zlim([-0.8,0]);

