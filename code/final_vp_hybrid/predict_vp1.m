function [ pf,R ] = predict_vp1( pf,u,dt,R)
global Param;

    uu=mvnrnd(u,Param.Qu); %random u
    R=get_R(pf.pose,uu,dt,R);
    %%%% --------------- calculate new pose -------------- %%%%
    x=pf.pose(1);
    y=pf.pose(2);
    phi=pf.pose(3);
    new_x=x+dt*(uu(1)*cos(phi)-(uu(1)/Param.L)*tan(uu(2))*(Param.a*sin(phi)+Param.b*cos(phi)));
    new_y=y+dt*(uu(1)*sin(phi)+(uu(1)/Param.L)*tan(uu(2))*(Param.a*cos(phi)-Param.b*sin(phi)));
    new_phi=minimizedAngle(phi+dt*(uu(1)/Param.L)*tan(uu(2)));
    %%%% ------------ add noise ------------------
    noise=mvnrnd(zeros(3,1),Param.Qf);
    new_x=new_x+noise(1);
    new_y=new_y+noise(2);
    new_phi=new_phi+noise(3);
    pf.pose(1)=new_x;
    pf.pose(2)=new_y;
    pf.pose(3)=new_phi;
    
              
end


function R=get_R(x,u,dt,R)    %called by get_Pi
global Param;
L=Param.L;
a=Param.a;
b=Param.b;
V=[dt*(cos(x(3))-1/L*tan(u(2))*(a*sin(x(3))+b*cos(x(3)))) -dt*u(1)/L/cos(u(2))^2*(a*sin(x(3))+b*cos(x(3)));
   dt*(sin(x(3))+1/L*tan(u(2))*(a*cos(x(3))-b*sin(x(3)))) dt*u(1)/L/cos(u(2))^2*(a*cos(x(3))-b*sin(x(3)));
   dt/L*tan(u(2))                                            dt*u(1)/L/cos(u(2))^2];
G=eye(3);
G(1,3)=dt*(-u(1)*sin(x(3))-u(1)/L*tan(u(2))*(a*cos(x(3))-b*sin(x(3))));
G(2,3)=dt*(u(1)*cos(x(3))+u(1)/L*tan(u(2))*(-a*sin(x(3))-b*cos(x(3))));

R=G*R*G'+V*Param.Qu*V'+Param.Qf;
end  %end for function get_R

