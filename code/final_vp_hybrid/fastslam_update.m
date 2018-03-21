function weight=fastslam_update(PI,ut,zt,Xtj,Hmj,Hxj,Qj,zj_bar,Rt,p)
%%%%
% the input are PI: line14, correspondence likelihood of all exsisting landmarks
%               ut: the current control 
%               zt: the current onservation
%               Xtj: line 12, the smaple pose of all exsisting landmarks, a vector
%               Hmj: line 8, should have the dimension as  2*(2*pf(p).Mn)
%               Hxj: line 7, should have the dimension as 2*(3*pf(p).Mn)
%               Qj: line 9, meansurement inforemation to jth landmark, should have the dimension 2*2*pf(p).Mn
%               zj_bar: line 13, meansurement prediction of jt landmark
%               Rt: the motion noise projected to state, should have the
%               dimension as 2*(2*pf(p).Mn)
%
%
%               Notice: everytime you call this function, it will change
%               the state of global pf.
global pf
global Param
%
% pf: particle 
%
p0=0.0005;
 Mt=[Param.alphas(1)*ut(1)^2+Param.alphas(2)*ut(2)^2 0 0;
            0 Param.alphas(3)*ut(2)^2+Param.alphas(4)*ut(1)^2+Param.alphas(4)*ut(3)^2 0;
            0 0 Param.alphas(1)*ut(3)^2+Param.alphas(2)*ut(2)^2];
PI(pf(p).Mn+1)=p0;% line 16, the likelyhood of new features, po is a parameter
[maxPI,c_bar]=max(PI); % line 17, C_bar: the most likelyhood landmarks, PI: the weight of pth particles to Nt-1 landmarks and new features
%
A=pf(p).Mn;
pf(p).Mn=max(pf(p).Mn,c_bar); % line 18: new number of landmarks
%
for j=1:pf(p).Mn
    if j==c_bar && c_bar==A+1
        ut=ut+sqrt([Mt(1,1);Mt(2,2);Mt(3,3)]).*randn(3,1); % line 21, sample u
        pf(p).pose=prediction(pf(p).pose,ut); % line 21, sample pose
        pf(p).lm(:,j)=obs_inv(zt,pf(p).pose); % line 22, initialize mean
        H=obs_J_map(pf(p).lm(:,j),pf(p).pose); % line 23
        
         Hm=get_Hm(zt,pf(p).pose);
        %pf(p).cov(:,2*j-1:2*j)=Hm*Param.R*Hm';% line 24, initialize covariance
        pf(p).cov(:,2*j-1:2*j)=inv(H)*Param.R*inv(H)';% line 24, initialize covariance
        pf(p).i(j)=1;
        weight=p0;% line 26, the output of this function.
    elseif j==c_bar && c_bar<=A
        pf(p).pose=Xtj(1:3,j);% line 28
        [zbar Hx Hm]=h_and_Jacob(pf(p).lm(:,j),pf(p).pose);
        K=pf(p).cov(:,2*j-1:2*j)*Hm'*Qj(1:2,2*j-1:2*j)^(-1);  % line 29
        diff=zt-zj_bar(1:2,j);
        diff(2)=minimizedAngle(diff(2));
        pf(p).lm(:,j)=pf(p).lm(:,j)+K*diff;  % line 30
        L=Hx*Rt*Hx'+Hm*pf(p).cov(:,2*j-1:2*j)*Hm'+Param.R;% line 33
        pf(p).cov(:,2*j-1:2*j)=(eye(2,2)-K*Hm)*pf(p).cov(:,2*j-1:2*j); % line 31
        pf(p).i(j)=pf(p).i(j)+1;
        diff=zt-zj_bar(1:2,j);
        diff(2)=minimizedAngle(diff(2));
        weight=det(2*pi*L)^(-0.5)*exp(-0.5*(diff)'*L^(-1)*(diff));% line 34, the output of this function.
    else
        pf(p).lm(:,j)=pf(p).lm(:,j); %line 36
        pf(p).cov(:,2*j-1:2*j)=pf(p).cov(:,2*j-1:2*j);    % line 37
    end
end


%%%%%
function Hm=get_Hm(z,x)
Hm=[cos(x(3)+z(2)) -z(1)*sin(x(3)+z(2));
    sin(x(3)+z(2)) z(1)*cos(x(3)+z(2))];



function state=prediction(state,motion)
state(3)=state(3)+motion(1);
state(1)=state(1)+motion(2)*cos(state(3));
state(2)=state(2)+motion(2)*sin(state(3));
state(3)=state(3)+motion(3);
state(3)=minimizedAngle(state(3));

%%%%%

function obs = observation(landmark,state)

dx = landmark(1) - state(1);
dy = landmark(2) - state(2);
dist = sqrt(dx^2 + dy^2);
obs = [dist ; minimizedAngle(atan2(dy, dx) - state(3))];


%%%%%

function H=obs_J_map(landmark,state)
dx = landmark(1) - state(1);
dy = landmark(2) - state(2);
dist = sqrt(dx^2 + dy^2);
H=[dx/dist, dy/dist;
    -dy/dist^2, dx/dist^2];

%%%%%

function landmark=obs_inv(z,state)

landmark(1)=state(1)+z(1)*cos(z(2)+state(3));
landmark(2)=state(2)+z(1)*sin(z(2)+state(3));

%%%%%


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