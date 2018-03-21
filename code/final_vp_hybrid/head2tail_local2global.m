function [mean, cov]=head2tail_local2global(l_mu,l_cov)   %r_mu and r_cov refer to 
%the robot pose and covariance in global coordinates immediately before
%submap initialization,  3*1, 3*3
%l_mu and l_cov refer to the landmark mean and covariance represented in
%the local coordinate   2n*1, 2n*2n
%mean and cov return the mean and cov of the landmarks in the global
%cooridinate   2n*1, 2n*2n
global State
x=State.Ekf.mu(1);
y=State.Ekf.mu(2);
t=State.Ekf.mu(3);   %t means theta

n=(length(l_mu)-3)/2;
m=(length(State.Ekf.mu)-3)/2;
% prepare the mean
mean=zeros(3+2*n+3+2*m,1);
% prepare the covariance
cov=zeros(3+2*n+3+2*m,3+2*n+3+2*m);
%
mean(1:3+2*m)=State.Ekf.mu;
%
%
%note: totally we have sixteen jacob
jacob_ogr_ogr=eye(3,3); % the Jacobian of old global robot to old global robot
jacob_ogl_ogl=eye(2*m,2*m);% the Jacobian of old global landmarks to old global landmarks
jacob_ogr_ogl=zeros(3,2*m);% the Jacobian of old global robot to old global landmarks
jacob_ogl_ogr=zeros(2*m,3);% the Jacobian of old global landmarks to old global robot
%
jacob_ogr_lr=zeros(3,3); % the Jacobian of old global robot to local robot
jacob_ogr_ll=zeros(3,2*n);% the Jacobian of old global robot to local landmark
%
jacob_ogl_lr=zeros(2*m,3); % the Jacobian of old global landmark to local robot
jacob_ogl_ll=zeros(2*m,2*n);% the Jacobian of old global landmark to local landmark
%
jacob_ngr_ogl=zeros(3,2*m);% the Jacobian of new global robot to old global landmark
jacob_ngl_ogl=zeros(2*n,2*m);% the Jacobian of new global landmark to old global landmark
%
jacob_ngr_ogr=[1,0,0;
               0,1,0;
               0,0,1];% the Jacobian of new global robot to old global robot
jacob_ngr_lr=[cos(State.Ekf.mu(3)),-sin(State.Ekf.mu(3)),0;
              sin(State.Ekf.mu(3)),cos(State.Ekf.mu(3)),0;
              0,0,1];% the Jacobian of new global robot to local robot
jacob_ngr_ll=zeros(3,2*n);% the Jacobian of new global robot to local landmark
%
jacob_ngl_lr=zeros(2*n,3);% the Jacobian of new global landmark to local robot
%
mean(3+2*m+1)=x+l_mu(1)*cos(t)-l_mu(2)*sin(t);
mean(3+2*m+2)=y+l_mu(1)*sin(t)+l_mu(2)*cos(t);
mean(3+2*m+3)=minimizedAngle(t+l_mu(3));
for i=1:n
    mean(2*i+3+2*m-1+3,1)=x+l_mu(2*i-1+3)*cos(t)-l_mu(2*i+3)*sin(t);
    mean(2*i+3+2*m+3,1)=y+l_mu(2*i-1+3)*sin(t)+l_mu(2*i+3)*cos(t);
end
%
for i=1:n
    
    jacob_ngl_ogr(2*i-1:2*i,1:3)=[1,0,-l_mu(2*i-1)*sin(t)-l_mu(2*i)*cos(t);
                            0,1, l_mu(2*i-1)*cos(t)-l_mu(2*i)*sin(t)];   % the Jacobian of new global landmark to old global robot
    
    jacob_ngl_ll(2*i-1:2*i,2*i-1:2*i)=[cos(t),-sin(t);
                                       sin(t),cos(t)];% the Jacobian of new global landmark to local landmark
    
end

J_total=zeros(3+2*m+3+2*n,3+2*m+3+2*n);
%
J_total(1:3,1:3)=jacob_ogr_ogr;
%
J_total(1:3,3+1:3+2*m)=jacob_ogr_ogl;
%
J_total(1:3,3+2*m+1:3+2*m+3)=jacob_ogr_lr;
%
J_total(1:3,3+2*m+3+1:3+2*m+3+2*n)=jacob_ogr_ll;
%%%%%
J_total(3+1:3+2*m,1:3)=jacob_ogl_ogr;
%
J_total(3+1:3+2*m,3+1:3+2*m)=jacob_ogl_ogl;
%
J_total(3+1:3+2*m,3+2*m+1:3+2*m+3)=jacob_ogl_lr;
%
J_total(3+1:3+2*m,3+2*m+3+1:3+2*m+3+2*n)=jacob_ogl_ll;
%%%%%
J_total(3+2*m+1:3+2*m+3,1:3)=jacob_ngr_ogr;
%
J_total(3+2*m+1:3+2*m+3,3+1:3+2*m)=jacob_ngr_ogl;
%
J_total(3+2*m+1:3+2*m+3,3+2*m+1:3+2*m+3)=jacob_ngr_lr;
%
J_total(3+2*m+1:3+2*m+3,3+2*m+3+1:3+2*m+3+2*n)=jacob_ngr_ll;
%%%%%
J_total(3+2*m+3+1:3+2*m+3+2*n,1:3)=jacob_ngl_ogr;
%
J_total(3+2*m+3+1:3+2*m+3+2*n,3+1:3+2*m)=jacob_ngl_ogl;
%
J_total(3+2*m+3+1:3+2*m+3+2*n,3+2*m+1:3+2*m+3)=jacob_ngl_lr;
%
J_total(3+2*m+3+1:3+2*m+3+2*n,3+2*m+3+1:3+2*m+3+2*n)=jacob_ngl_ll;
%


cov(1:3,1:3)=State.Ekf.Sigma(1:3,1:3);
cov(1:3,3+1:3+2*m)=State.Ekf.Sigma(1:3,3+1:3+2*m);
cov(1:3,3+2*m+1:3+3+2*m+2*n)=State.Ekf.Sigma(1:3,1:3)'*J_total(3+2*m+1:3+2*m+3+2*n,1:3)';
cov(3+1:3+2*m,1:3)=State.Ekf.Sigma(1:3,3+1:3+2*m)';
cov(3+1:3+2*m,3+1:3+2*m)=State.Ekf.Sigma(3+1:3+2*m,3+1:3+2*m);
cov(3+1:3+2*m,3+2*m+1:3+2*m+2*n+3)=State.Ekf.Sigma(1:3,3+1:3+2*m)'*J_total(3+2*m+1:3+2*m+3+2*n,1:3)';
cov(3+2*m+1:3+2*m+3+2*n,1:3)=cov(1:3,3+2*m+1:3+3+2*m+2*n)';
cov(3+2*m+1:3+2*m+3+2*n,3+1:3+2*m)=cov(3+1:3+2*m,3+2*m+1:3+2*m+2*n+3)';
cov(3+2*m+1:3+2*m+3+2*n,3+2*m+1:3+2*m+2*n+3)=J_total(3+2*m+1:3+2*m+3+2*n,1:3)*State.Ekf.Sigma(1:3,1:3)*J_total(3+2*m+1:3+2*m+3+2*n,1:3)'+...
    J_total(3+2*m+1:3+2*m+3+2*n,3+2*m+1:3+2*m+3+2*n)*l_cov*J_total(3+2*m+1:3+2*m+3+2*n,3+2*m+1:3+2*m+3+2*n)';

end

