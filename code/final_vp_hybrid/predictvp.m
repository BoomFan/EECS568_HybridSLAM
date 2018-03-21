function [ new,pf,chat,p0,zt_zhat,Rt ] = predictvp( pf,u,dt,zt,Rt)
global Param;
xhat=g(pf.pose,u,dt);
Rt=get_R(pf.pose,u,dt,Rt);
p0=1.5;
p00=0.003;
spur=0;
if pf.Mn==0
    new=1;
    uu=mvnrnd(u,Param.Qu); %random u
    pf.pose=g(pf.pose,uu,dt);
    pf.pose=pf.pose+mvnrnd([0;0;0],Param.Qf);
    
    Pi(pf.Mn+1)=p0;
%     pf.Mn=pf.Mn+1;
    zt_zhat=zeros(2,1);
    chat=1;
else
    
    
    for j=1:pf.Mn
        [zbar Hx Hm]=h_and_Jacob(pf.lm(:,j),xhat);   %line6,7,8
        Q=Param.R+Hm*pf.cov(:,2*j-1:2*j)*Hm';
        Q=Q*10;
       
        Cov_proposal=inv(Hx'*inv(Q)*Hx+Rt^-1);
        zt_zbar=zt-zbar;
        zt_zbar(2)=minimizedAngle(zt_zbar(2));   %minimize angle
        mu_proposal=Cov_proposal*Hx'*inv(Q)*(zt_zbar)+xhat';                %line 11
        mu_proposal(3)=minimizedAngle(mu_proposal(3));    %minimize angle
        xtj(1:3,j) =mu_proposal+chol(Cov_proposal)'*randn(3,1);

        zhat=h(pf.lm(:,j),xtj(1:3,j));
        zt_zhat(:,j)=zt-zhat;
        zt_zhat(2,j)=minimizedAngle(zt_zhat(2,j));
        Q=5*Q;
        Pi(j)=1/sqrt(det(2*pi*Q))*exp(-0.5*zt_zhat(:,j)'*inv(Q)*zt_zhat(:,j));
        
    end
    %add
%     test=max(Pi);
    Pi(pf.Mn+1)=p0;
    [maxPI,chat]=max(Pi);
    A=pf.Mn;
%     pf.Mn=max(pf.Mn,chat);
    if A==max(pf.Mn,chat)
        new=0;
        pf.pose=xtj(1:3,chat);
    else
        new=1;
        uu=mvnrnd(u,Param.Qu); %random u
        pf.pose=g(pf.pose,uu,dt)+mvnrnd([0;0;0],Param.Qf);
%         if max(Pi)>p00
%             spur=1;
%         end
    end

    
end



end








































function xhat = g(x,u,dt)    %called by get_Pi
global Param
L=Param.L;
a=Param.a;
b=Param.b;
xhat(1)=x(1)+dt*(u(1)*cos(x(3))-u(1)/L*tan(u(2))*(a*sin(x(3))+b*cos(x(3))));
xhat(2)=x(2)+dt*(u(1)*sin(x(3))+u(1)/L*tan(u(2))*(a*cos(x(3))-b*sin(x(3))));
xhat(3)=minimizedAngle(x(3)+dt*u(1)/L*tan(u(2)));
end  %end for function_g



function [zbar Hx Hm]=h_and_Jacob(mu,xhat)    %called by get_Pi
dx=mu(1)-xhat(1);
dy=mu(2)-xhat(2);
q=dx^2+dy^2;
zbar=[sqrt(q);
    atan2(dy,dx)-xhat(3)]; %prediction of observation of bearing
zbar(2)=minimizedAngle(zbar(2));

Hx=1/q*[-sqrt(q)*dx,-sqrt(q)*dy,0;
    dy,         -dx,      -q];

Hm=1/q*[sqrt(q)*dx,sqrt(q)*dy;
    -dy,      dx];
end  %end for function_h_and_Jacob

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

function zbar=h(mu,xhat)     %called by get_Pi
dx=mu(1)-xhat(1);
dy=mu(2)-xhat(2);
q=dx^2+dy^2;
zbar=[sqrt(q);
    atan2(dy,dx)-xhat(3)]; %prediction of observation of bearing
zbar(2)=minimizedAngle(zbar(2));
end  %end for function_h_and_Jacob

function state=prediction(state,motion)
state(3)=state(3)+motion(1);
state(1)=state(1)+motion(2)*cos(state(3));
state(2)=state(2)+motion(2)*sin(state(3));
state(3)=state(3)+motion(3);
state(3)=minimizedAngle(state(3));
end

