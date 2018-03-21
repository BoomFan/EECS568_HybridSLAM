function [ new,pf,chat,p0,zt_zhat,Rt ] = newpredict( pf,u,zt )
global Param;
xhat=g(pf.pose,u);
Rt=get_R(pf.pose,u);

if pf.Mn==0
    new=1;
    Mt=[Param.alphas(1)*u(1)^2+Param.alphas(2)*u(2)^2 0 0;
        0 Param.alphas(3)*u(2)^2+Param.alphas(4)*u(1)^2+Param.alphas(4)*u(3)^2 0;
        0 0 Param.alphas(1)*u(3)^2+Param.alphas(2)*u(2)^2];
    uu=mvnrnd(u,Mt); %random u
    pf.pose=g(pf.pose,uu);
    p0=0.0005;
    Pi(pf.Mn+1)=p0;
    pf.Mn=pf.Mn+1;
    zt_zhat=zeros(2,1);
    chat=1;
else
    
    
    for j=1:pf.Mn
        [zbar Hx Hm]=h_and_Jacob(pf.lm(:,j),xhat);   %line6,7,8
        Q=Param.R+Hm*pf.cov(:,2*j-1:2*j)*Hm';
        
        if abs(det(Rt))<=10e-9
            ut=u;
            Mt=[Param.alphas(1)*ut(1)^2+Param.alphas(2)*ut(2)^2 0 0;
            0 Param.alphas(3)*ut(2)^2+Param.alphas(4)*ut(1)^2+Param.alphas(4)*ut(3)^2 0;
            0 0 Param.alphas(1)*ut(3)^2+Param.alphas(2)*ut(2)^2];
        ut=ut+sqrt([Mt(1,1);Mt(2,2);Mt(3,3)]).*randn(3,1);
        
            xtj(1:3,j) =prediction(pf.pose,ut);
        else
            Cov_proposal=inv(Hx'*inv(Q)*Hx+Rt^-1);
            zt_zbar=zt-zbar;
            zt_zbar(2)=minimizedAngle(zt_zbar(2));   %minimize angle
            mu_proposal=Cov_proposal*Hx'*inv(Q)*(zt_zbar)+xhat;                %line 11
            mu_proposal(3)=minimizedAngle(mu_proposal(3));    %minimize angle
            xtj(1:3,j) =mu_proposal+chol(Cov_proposal)'*randn(3,1);
        end
        zhat=h(pf.lm(:,j),xtj(1:3,j));
        zt_zhat(:,j)=zt-zhat;
        zt_zhat(2,j)=minimizedAngle(zt_zhat(2,j));
        Pi(j)=1/sqrt(det(2*pi*Q))*exp(-0.5*zt_zhat(:,j)'*inv(Q)*zt_zhat(:,j));
        
    end
    p0=0.00005;
    Pi(pf.Mn+1)=p0;
    [maxPI,chat]=max(Pi);
    A=pf.Mn;
    pf.Mn=max(pf.Mn,chat);
    if A==pf.Mn
        new=0;
        pf.pose=xtj(1:3,chat);
    else
        new=1;
        Mt=[Param.alphas(1)*u(1)^2+Param.alphas(2)*u(2)^2 0 0;
            0 Param.alphas(3)*u(2)^2+Param.alphas(4)*u(1)^2+Param.alphas(4)*u(3)^2 0;
            0 0 Param.alphas(1)*u(3)^2+Param.alphas(2)*u(2)^2];
        uu=mvnrnd(u,Mt); %random u
        pf.pose=g(pf.pose,uu);
    end
    
end


end








































function xhat = g(x,u)    %called by get_Pi
xhat=x+[u(2)*cos(x(3)+u(1));
    u(2)*sin(x(3)+u(1));
    u(1)+u(3)];
xhat(3)=minimizedAngle(xhat(3));
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

function R=get_R(x,u)    %called by get_Pi
global Param;
Vt=[-u(2)*sin(x(3)+u(1)), cos(x(3)+u(1)), 0;
    u(2)*cos(x(3)+u(1)) sin(x(3)+u(1)) 0;
    1 0 1];   %paetial xt+1 w.r.t. ut
Mt=[Param.alphas(1)*u(1)^2+Param.alphas(2)*u(2)^2 0 0;
    0 Param.alphas(3)*u(2)^2+Param.alphas(4)*u(1)^2+Param.alphas(4)*u(3)^2 0;
    0 0 Param.alphas(1)*u(3)^2+Param.alphas(2)*u(2)^2];

R=Vt*Mt*Vt.';
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

