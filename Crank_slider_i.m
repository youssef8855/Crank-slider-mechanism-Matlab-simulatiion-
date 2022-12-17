% Crank slider mechanism
clear all
close all
clc

l2=0.2; l3=0.3; omg2=4; tho2=14.48*pi/180; h=.05; 
dt=0.1;% Step size
t_end=4; t_start=0; 
n_sol=(t_end-t_start)/dt+1; %Number of Steps
error_tol=1.0E-6;

Res_mat=zeros(n_sol,37);

q_num=[0 0 0 0.5*l2*cos(tho2) .5*h tho2 l2*cos(tho2)+0.5*l3 h 0 l2*cos(tho2)+l3 h 0];
qd_num=zeros(1,12);
qdd_num=zeros(1,12);

%% constraint 
syms RX1 RY1 TH1 RX2 RY2 TH2 RX3 RY3 TH3 RX4 RY4 TH4 t
syms RX1d RY1d TH1d RX2d RY2d TH2d RX3d RY3d TH3d RX4d RY4d TH4d 
q=[RX1 RY1 TH1 RX2 RY2 TH2 RX3 RY3 TH3 RX4 RY4 TH4];
qd=[RX1d RY1d TH1d RX2d RY2d TH2d RX3d RY3d TH3d RX4d RY4d TH4d ];

C=[ RX1;
    RY1;
    TH1;
    RX2-0.5*l2*cos(TH2);
    RY2-0.5*l2*sin(TH2);
    RX2+0.5*l2*cos(TH2)-RX3+0.5*l3*cos(TH3);
    RY2+0.5*l2*sin(TH2)-RY3+0.5*l3*sin(TH3);
    RX3+0.5*l3*cos(TH3)-RX4;
    RY3+0.5*l3*sin(TH3)-RY4;
    RY4;
    TH4;
    TH2-omg2*t-tho2];
%% For position
Cq=jacobian(C,q);
% For velocity you need also
Ct=diff(C,t);
% For acceleration you need in addition to Cq the following
Ctt=diff(Ct,t);
Cqt=diff(Cq,t);

Cq_qd=Cq*qd.';
Cq_qdq=jacobian(Cq_qd,q);

Qd=-Cq_qdq*qd.'-2*Cqt*qd.'-Ctt;

%% At t=0
% For verification
C_num1=subs(C,q,q_num);
C_num2=subs(C_num1,t,0);
Cq_num=subs(Cq,q,q_num);
qd_num=-(Cq_num\Ct)';

Qd_num1=subs(Qd,q,q_num);
Qd_num2=subs(Qd_num1,qd,qd_num);
qdd_num=Cq_num\Qd_num2;

Res_mat(1,2:13)=q_num;
Res_mat(1,14:25)=qd_num;
Res_mat(1,26:37)=qdd_num;
Res_mat(1:n_sol,1)=t_start:dt:t_end;



for i_res=2:n_sol
    t_num=Res_mat(i_res,1);
    q_num_n=Res_mat(i_res-1,2:13)+dt*Res_mat(i_res-1,14:25);
    error1=1.0;
    while abs(error1)>error_tol,   
    C_num1=subs(C,q,q_num_n);
    C_num2=subs(C_num1,t,t_num);
    Cq_num=subs(Cq,q,q_num_n);
    C_num2 = vpa(C_num2);
    C_num2 = simplify(C_num2);
    Cq_num = vpa(Cq_num);
    Cq_num = simplify(Cq_num);
     
    d_q_num_n=-(Cq_num\C_num2)';
    q_num_np1=q_num_n+d_q_num_n;
    error1=eval(norm(C_num2));
    error2=eval(norm(d_q_num_n));
    q_num_n=q_num_np1;
    
    end
    
    Cq_num=subs(Cq,q,q_num_n);
    qd_num=-(Cq_num\Ct)';

    Qd_num1=subs(Qd,q,q_num_n);
    Qd_num2=subs(Qd_num1,qd,qd_num);
    qdd_num=Cq_num\Qd_num2;
    
    Res_mat(i_res,2:13)=q_num_n;
    Res_mat(i_res,14:25)=qd_num;
    Res_mat(i_res,26:37)=qdd_num;
    
end
%%
%required(1)
figure(1);
plot(Res_mat(:,1), Res_mat(:,11));
title('(1)Horizontal position of slider');
xlabel('t');
ylabel('x');

%required(2)
figure(2);
plot(Res_mat(:,1),Res_mat(:,23));
title('(2)Horizontal velocity of slider');
xlabel('t');
ylabel('v');

%required(3)
figure(3);
Point_A= sin(Res_mat(:,7))*l2;
plot(Res_mat(:,1),Point_A);
title('(3)Vertical position of point A');
xlabel('t');
ylabel('yA');
%required(4)
PX= Res_mat(:,8)+ 0.1*cos(Res_mat(:,10))-0.2*sin(Res_mat (:,10));
PY= Res_mat(:,9)+0.1*sin(Res_mat(:,10))+0.2*cos(Res_mat(:,10));
figure(4);
plot(PX,PY);
title(' (4)Trace of  point (0.1, 0.2) defined in body 3 coordinate system.');
xlabel('X ');
ylabel('Y ');




