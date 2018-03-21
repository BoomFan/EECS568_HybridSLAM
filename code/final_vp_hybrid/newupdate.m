function [ w,pf ] = newupdate( new,pf,chat,p0,z,zt_zhat,Rt )
global Param;
j=chat;
if new
    if pf.Mn~=chat  % a test
        'error'
    end
    pf.lm(:,j)=obs_inv(z,pf.pose);
    [dummy1 dummy2 Hm]=h_and_Jacob(pf.lm(:,j),pf.pose);
    pf.cov(:,2*j-1:2*j)=inv(Hm)*Param.R*inv(Hm)';
    pf.i(j)=1;
    w=p0;

else
    [zbar Hx Hm]=h_and_Jacob(pf.lm(:,j),pf.pose);
    Q=Param.R+Hm*pf.cov(:,2*j-1:2*j)*Hm';
    K=pf.cov(:,2*j-1:2*j)*Hm'*Q^-1;  %line 29;
    
    A=pf.cov(:,2*j-1:2*j); % for line 33 use
    
    pf.lm(:,j)=pf.lm(:,j)+K*zt_zhat(:,j);  % line 30
    pf.cov(:,2*j-1:2*j)=(eye(2,2)-K*Hm)*pf.cov(:,2*j-1:2*j);
    L=Hx*Rt*Hx'+Hm*A*Hm'+Param.R;
    w=det(2*pi*L)^(-0.5)*exp(-0.5*zt_zhat(:,j)'*L^(-1)*zt_zhat(:,j));
    pf.i(j)=pf.i(j)+1;
end
end






function landmark=obs_inv(z,state)

landmark(1)=state(1)+z(1)*cos(z(2)+state(3));
landmark(2)=state(2)+z(1)*sin(z(2)+state(3));
end


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