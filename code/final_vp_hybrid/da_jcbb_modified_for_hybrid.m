function Li=da_jcbb_modified_for_hybrid(z, cov_lm)   %cov_lm is the covariance matrix of the new landmarks   
% perform joint-compatability branch and bound data association

global Param;
global State;
global Best_pair;   %global variable shared with subfunction JCBB
global z_JCBB;  %global variable shared with subfunction JCBB
global P;
z_JCBB=z;
P=cov_lm;
m=length(z)/2;   %length of observation
Best_pair=zeros(1,m);  %initialize best as m null hypothesis
global NL
NL=State.Ekf.nL;  
if(NL==0)  %the first time step,all observations are used for state augmentation
    for i=1:m
        Li(3,i)=0;  %indicate new landmark
    end
    
else  %else for if(State.Ekf.nL==0)  not the first time step
    
    JCBB([],1);    %find pairings for feature Ei, i<=m, use null hypothesis and the first observation to initialize
    %after calling recursive function JCBB, Best_pair should have some nonzero hypothesis
    %initilize the zero in Best_pair here
    for i=1:m
        if Best_pair(i)==0
            
            Li(3,i)=0;  %new landmark
            
            for j=1:State.Ekf.nL
                if unitary_compatible(i,j,0.99999)
                    Best_pair(i)=-1;  %indicate this observation is individually compatible but jointly incompatible, can be spurious
                    %it should be thrown
                    Li(3,i)=-1;
                    break
                end
            end


        else
        Li(1,i)=3+2*Best_pair(i)-1;
        Li(2,i)=3+2*Best_pair(i);
        Li(3,i)=1;  %known landmark
        end
    end
    
end % end for if (State.Ekf.nL==0)

end %end for da_jcbb





function JCBB(H,i)
global Param;
global State;
global Best_pair;
global z_JCBB;
global P;
m=length(z_JCBB)/2;
if (i>m)
    
    if sum(H>0)>sum(Best_pair>0)
    Best_pair=H;
    end  

else % else for if (i>m)
    
    for j=1:State.Ekf.nL
        
        if unitary_compatible(i,j,0.999)

        if joint_compatible(H,i,j)
      
    JCBB([H j],i+1);
           
        end
        end
    end  %end for
    
    if sum(H>0)+m-i>sum(Best_pair>0)
        JCBB([H 0],i+1);
    end
end  %end for if (i>m)

end  %end JCBB

function y=unitary_compatible(i,j,p0)
global Param;
global State;
global z_JCBB;
global P;
D=Mahalanobis2(z_JCBB(2*i-1:2*i),i,j);
if D<=chi2inv(p0,2)
    y=1;
else
    y=0;
end
end


function y=joint_compatible(H,i,j)   %j is the indices in ekf, i is the indices of z_jcbb
global Param;
global State;
global z_JCBB;
global P;
N=length(H);
m=length(z_JCBB)/2;
NN=sum(H>0);
       HH=zeros(2*NN+2,2*(State.Ekf.nL+m));  %all the non-null hypothesis and i,j
       M=1;
       
        Hlow=[1 0 , -1 0;
       0 1 , 0  -1];
       for k=1:N
           if (H(k)~=0)
        HH(M:M+1,2*H(k)-1:2*H(k))=Hlow(1:2,1:2);
        HH(M:M+1,2*State.Ekf.nL+2*i-1:2*State.Ekf.nL+2*i)=Hlow(1:2,3:4);
%         hH(M:M+1,1)=zt;
         z_diff(M:M+1,1)=[State.Ekf.mu(3+2*H(k)-1)-z_JCBB(2*k-1);
                          State.Ekf.mu(3+2*H(k))-z_JCBB(2*k)];
        M=M+2;
       end  %end if (H(k)!=0)
       end %end for k=1:N
       if NN==0
z_diff=[
    State.Ekf.mu(3+2*j-1)-z_JCBB(2*i-1);
    State.Ekf.mu(3+2*j)-z_JCBB(2*i)];    %final pairing i,j
    HH(1:2,2*j-1:2*j)=Hlow(1:2,1:2);
    HH(1:2,2*State.Ekf.nL+2*i-1:2*State.Ekf.nL+2*i)=Hlow(1:2,3:4);
       else
         z_diff=[z_diff;
    State.Ekf.mu(3+2*j-1)-z_JCBB(2*i-1);
    State.Ekf.mu(3+2*j)-z_JCBB(2*i)];    %final pairing i,j

    
    HH(end-1:end,2*j-1:2*j)=Hlow(1:2,1:2);
    HH(end-1:end,2*State.Ekf.nL+2*i-1:2*State.Ekf.nL+2*i)=Hlow(1:2,3:4);  
       end
        PCH=zeros(2*(State.Ekf.nL+m),2*(State.Ekf.nL+m));
        PCH(1:2*State.Ekf.nL,1:2*State.Ekf.nL)=State.Ekf.Sigma(4:end,4:end);
        PCH(2*State.Ekf.nL+1:end,2*State.Ekf.nL+1:end)=P;
        CH=HH*PCH*HH.';
if det(CH)<=10^(-10)
    y=0;
else
        DH2=z_diff.'*inv(CH)*z_diff;
        
                        if DH2<0|imag(DH2)~=0
            'Warning:M-distance(DH2) negative or imaginary'
                        end

                
        if DH2<=chi2inv(0.999,2*NN+2)  %here N is the total length of H but NN is the # of null hypothesis
            y=1;
        else
            y=0;
        end  %end for if
end
        
end  %end for the function
        
        
        
        
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