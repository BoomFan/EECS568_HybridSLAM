function [Pi,xtj,Hx_final,Hm_final,Qj_final,zj_bar,R]=get_predict(pf,u,zt)  %here pf is a local quantity, which is pf(k) in the 
    global Param;
    xhat=g(pf.pose,u);               %line5, xhat does not depend on j
    if pf.Mn==0
        
        
        
    else
    for j=1:pf.Mn

    [zbar Hx Hm]=h_and_Jacob(pf.lm(:,j),xhat);      %line6,7,8
    Hx_final(1:2,3*j-2:3*j)=Hx;
    Hm_final(1:2,2*j-1:2*j)=Hm;
    Q=Param.R+Hm*pf.cov(:,2*j-1:2*j)*Hm.';    %line 9,Param.R is Qt in the book
    Qj_final(1:2,2*j-1:2*j)=Q;
     R=get_R(pf.pose,u);
    if abs(det(R))<=10^-7
        ut=u;
         Mt=[Param.alphas(1)*ut(1)^2+Param.alphas(2)*ut(2)^2 0 0;
            0 Param.alphas(3)*ut(2)^2+Param.alphas(4)*ut(1)^2+Param.alphas(4)*ut(3)^2 0;
            0 0 Param.alphas(1)*ut(3)^2+Param.alphas(2)*ut(2)^2];
        ut=ut+sqrt([Mt(1,1);Mt(2,2);Mt(3,3)]).*randn(3,1);
    xtj(1:3,j) =prediction(pf.pose,ut);
    else
    Cov_proposal=inv(Hx.'*inv(Q)*Hx+R^-1);   %line10
    zt_zbar=zt-zbar;
    zt_zbar(2)=minimizedAngle(zt_zbar(2));   %minimize angle
    mu_proposal=Cov_proposal*Hx.'*inv(Q)*(zt_zbar)+xhat;                %line 11
    mu_proposal(3)=minimizedAngle(mu_proposal(3));    %minimize angle
    L = chol(Cov_proposal);
    L=L.';
    xtj(1:3,j) = mu_proposal + L*randn(3,1);   %sample from multivariate Gaussian,need to be checked ,line 12
    end

    zhat=h(pf.lm(:,j),xhat);  %line13, is is xhat or pf.pose? Need to be checked 
    zj_bar(1:2,j)=zhat;
    zt_zhat=zt-zhat;
    zt_zhat(2)=minimizedAngle(zt_zhat(2));
    Pi(j)=1/sqrt(det(2*pi*Q))*exp(-0.5*zt_zhat.'*inv(Q)*zt_zhat);   %line14
    end   %end for j=1:pf(k).Mn  line15
    end
    end  %end for function get_Pi
    
  
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
        
        R=Vt*Mt*Vt.';   %projection of noise parameters from control space to robot pose space
%         if det(R)==0
%         R(1,2)=0;
%         R(2,1)=0;
%         else
%         end
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