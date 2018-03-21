function Li=da_dgnn_modified_for_hybrid(z, cov_lm)   %cov_lm is the covariance matrix of the new landmarks   
% perform joint-compatability branch and bound data association

global Param;
global State;
global z_dgnn;  %global variable shared with subfunction JCBB
global P;
P=cov_lm;
m=length(z)/2;   %length of observation
z_dgnn=z;
global NL
NL=State.Ekf.nL;  
if(NL==0)  %the first time step,all observations are used for state augmentation
    for i=1:m
        Li(3,i)=0;  %indicate new landmark
    end
    
    else  %else for if(State.Ekf.nL==0)  not the first time step
    
    for i=1:m
            nearest=Mahalanobis2(z_dgnn(2*i-1:2*i),i,1);
            nearest_indice=1;
            for j=2:State.Ekf.nL
                temp=Mahalanobis2(z_dgnn(2*i-1:2*i),i,j);
                if temp<nearest
                    nearest=temp;
                    nearest_indice=j;
                end
            end
            if nearest<chi2inv(0.998,2)
                    Li(1,i)=3+2*nearest_indice-1;  
                    Li(2,i)=3+2*nearest_indice;%indicate this observation is individually compatible but jointly incompatible, can be spurious
                    %it should be thrown
                    Li(3,i)=1;
            end
             if nearest>chi2inv(0.9999999999,2)
                    Li(3,i)=0;
                end
                if nearest<=chi2inv(0.999999,2)&nearest>=chi2inv(0.998,2)
                   Li(3,i)=-1; 
            end
                



%         else
%         Li(1,i)=3+2*Best_pair(i)-1;
%         Li(2,i)=3+2*Best_pair(i);
%         Li(3,i)=1;  %known landmark
%         end
    end
    
end % end for if (State.Ekf.nL==0)

end %end for da_jcbb



function y=Mahalanobis2(z,i,k)   %need to be modified,z is the i_th observation mean, i is the indice of the observation, k is the indice of the old landmark
global Param;
global State;
global P;
%         dx=State.Ekf.mu(3+2*k-1)-State.Ekf.mu(1);
%         dy=State.Ekf.mu(3+2*k)-State.Ekf.mu(2);
%         q=dx^2+dy^2;
              
%         switch lower(State.Ekf.simchoice)
%     case 'sim'
%  zt=[sqrt(q);
%      atan2(dy,dx)-State.Ekf.mu(3)]; %prediction of observation of bearing
%         zt(2)=minimizedAngle(zt(2));
%             case 'vp'
%                  zt=[sqrt(q);
%      atan2(dy,dx)-State.Ekf.mu(3)+pi/2]; %prediction of observation of bearing
%         zt(2)=minimizedAngle(zt(2));
%         end
        

%         Hlow=1/q*[-sqrt(q)*dx,-sqrt(q)*dy,0,sqrt(q)*dx,sqrt(q)*dy;
%           dy,         -dx,        -q,-dy,      dx];
Hlow=[1 0 , -1 0;
      0  1,  0 -1];
        PBF(1:2,1:2)=State.Ekf.Sigma(3+2*k-1:3+2*k,3+2*k-1:3+2*k);
        PBF(1:2,3:4)=zeros(2,2);
        PBF(3:4,1:2)=zeros(2,2);   %because local map is assumed to be uncorrelated with the global map
        PBF(3:4,3:4)=P(2*i-1:2*i,2*i-1:2*i);
        PBF=Hlow*PBF*Hlow.';
       % zmin=z-zt;
       % zmin(2)=minimizedAngle(zmin(2));
       z_diff=[State.Ekf.mu(3+2*k-1)-z(1);
           State.Ekf.mu(3+2*k)-z(2)];
        D2=z_diff.'*inv(PBF)*z_diff;
                if D2<0|imag(D2)~=0
            'Warning:M-distance negative or imaginary'
                end
        y=D2;
end