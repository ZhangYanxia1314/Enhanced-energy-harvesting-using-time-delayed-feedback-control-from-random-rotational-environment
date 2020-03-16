%采用迭代法计算能量依赖的频率w(H)
clc;
clear all;
syms newW;
global H w delta1 delta3 k alpha
delta1=3;
delta3=3;
k=0.3;
alpha=0.05;
hx=81;
hy=81;
hd=0.049;
x=-2+(0:hx-1)*hd;
y=-2+(0:hy-1)*hd;
n=length(y);
Wi=zeros(n,n);
for i=1:n
    for j=1:n
        tic
        Itn=1;w=0.001;
        while 1
            H=Hh(w,x(i),y(j));
            NH(Itn)=H;
            if H>=0
                Th=4*tt(0,xb(H,w));
                newW=2*pi/Th;
            else
                Th=2*tt(xa(H,w),xb(H,w));
                newW=2*pi/Th;
            end
            wt(Itn)=w;
            newWt(Itn)=newW;
            if Itn>500 %迭代次数
                m=abs(wt-newWt);
                ss=find(m==min(m));
                Wi(i,j)=newWt(ss);
                NHi(i,j)=NH(ss);               
                break;
            end
            w=w+0.01;
            Itn=Itn+1;
        end
        [n,i,j,Itn,w,newW,Wi(i,j)]
        toc
    end
end

Wisb=real(Wi);
figure;
subplot(1,2,1),mesh(x,y,NHi')
xlabel('$x$','Interpreter','latex');
ylabel('$\dot x$','Interpreter','latex');
zlabel('$H$','Interpreter','latex');
subplot(1,2,2),mesh(x,y,Wisb');
title({['\delta_1=',num2str(delta1),' \delta_3=',num2str(delta3),' \kappa=',num2str(k),' \alpha=',num2str(alpha)]})
xlabel('$x$','Interpreter','latex');
ylabel('$\dot x$','Interpreter','latex');
zlabel('$\omega(H)$','Interpreter','latex');

filename1=strcat('H_delta1=',num2str(delta1),'_delta3=',num2str(delta3),'_k=',num2str(k),'_alpha=',num2str(alpha),'.mat');
save(filename1,'NHi');
filename2=strcat('Wi_delta1=',num2str(delta1),'_delta3=',num2str(delta3),'_k=',num2str(k),'_alpha=',num2str(alpha),'.mat');
save(filename2,'Wisb');

%%%平滑频率Fig图底部毛刺高误差数值，优化频率结果
Win=Wisb;
for i=1:n
    for j=1:n
        if Win(i,j)<0.7 && abs(x(i))<1
            Win(i,j)=0.7;
        elseif Win(i,j)<0.7 && abs(x(i))>1
            Win(i,j)=0.75;
        end
    end
end
figure;
mesh(x,y,Win')
xlabel('$x$','Interpreter','latex');
ylabel('$\dot x$','Interpreter','latex');
zlabel('$\omega(H)$','Interpreter','latex');

filename1=strcat('WH_delta1=',num2str(delta1),'_delta3=',num2str(delta3),'_k=',num2str(k),'_alpha=',num2str(alpha),'.mat');
save(filename1,'Win');
