function Ekf_Hybrid_update(mu,cov,Li)    %mu and cov are the mean and covariance of the features in the local map represented in the global coordinates
%the third row of Li is the coorespondence associating each landmark in mu to the features
%in the old global map
%Li(3,i)=0 means new landmark, should be added directly to the global map
%Li(3,i)=1 means old landmark, should be fused to the old global landmark through modified Ekf update
%to produce a single consistent landmark
%Li(3,i)=-1 means JCBB regards this landmark as a spurious one, should
%throw it away
%mu is (3+2*m)*1,  the first three is the robot pose, m is the # of the
%landmarks in local map
global State;
mu(3)=minimizedAngle(mu(3));


m=(length(mu)-3)/2;   %total # of the newly observed landmarks
%First initialize new landmark


nL=State.Ekf.nL;
%Extract_indice=[4+2*nL:6+2*nL,4:3+2*nL];   %first include the original landmark
Extract_indice=2*nL+1:2*nL+3;   %first include the original pose and landmarks
Extract_indice=[Extract_indice, 1:2*nL];
N=sum(Li(3,:)>0);   %gives the # of landmarks that needs to be paired 

P=cov(4:end,4:end);  
xhat=[State.Ekf.mu(4:end);
            mu];
% if nL==0
%     State.Ekf.mu=mu;
%     State.Ekf.Sigma=P;
%     State.Ekf.nL=m;
% end


if N~=0

        

%if N~=0    %need to do update
HH=zeros(2*N,length(State.Ekf.mu)+length(mu)-3);
M=0;

for k=1:m
    if Li(3,k)==1   %indicate the k th landmark in mu need to be fused with the global landmark
        M=M+1;
    HH(2*M-1:2*M,Li(1,k)-3:Li(2,k)-3)=eye(2);  %revised@4_10
    HH(2*M-1:2*M,2*nL+3+2*k-1:2*nL+3+2*k)=-eye(2);
    hH(2*M-1,1)=State.Ekf.mu(Li(1,k))-mu(3+2*k-1);
    hH(2*M,1)=State.Ekf.mu(Li(2,k))-mu(3+2*k);
    end
end

if M~=N
    'Warning: # of landmark that needs fusing mismatch, ekf_hybrid_update line 34'
end


        %modified EKF update according to the paper
        K=P*HH.'*inv(HH*P*HH.');

        xhat=xhat-K*hH;
        xhat(2*nL+3)=minimizedAngle(xhat(2*nL+3));
%         K*hH
%         'update line 66'
%         xhat(3)=minimizedAngle(xhat(3));
%         xhat(length(State.Ekf.mu)+3)=minimizedAngle(xhat(length(State.Ekf.mu)+3));

        temp=K*HH;
        temp=eye(size(temp,1))-temp;
        P=temp*P;
        
%end %end if N~=0

for k=1:m
    if Li(3,k)==0   %new landmark
        State.Ekf.nL=State.Ekf.nL+1;   %update the # of landmarks
%          State.Ekf.mu(3+2*State.Ekf.nL-1)=mu(3+2*k-1);
%          State.Ekf.mu(3+2*State.Ekf.nL)=mu(3+2*k);
        Extract_indice=[Extract_indice,2*nL+3+2*k-1,2*nL+3+2*k];  %those need to be extracted from the fianl matrix
    end % end if Li(3,k)==0
end

% if nL~=0
%     xhat
% end
State.Ekf.Sigma=P([Extract_indice],[Extract_indice]); %another half has been discarded
State.Ekf.mu=xhat([Extract_indice]);   %another half has been discarded
else
    Extract_indice=2*nL+1:2*nL+3;
    Extract_indice=[Extract_indice, 1:2*nL];
    
for k=1:m
    if Li(3,k)==0   %new landmark
        State.Ekf.nL=State.Ekf.nL+1;   %update the # of landmarks
        Extract_indice=[Extract_indice,2*nL+3+2*k-1,2*nL+3+2*k];  %those need to be extracted from the fianl matrix
    end % end if Li(3,k)==0
end
State.Ekf.Sigma=P([Extract_indice],[Extract_indice]); %another half has been discarded
State.Ekf.mu=xhat([Extract_indice]);   %another half has been discarded
end  %end for if nL==0
if 3+2*State.Ekf.nL~=length(State.Ekf.mu)
    3+2*State.Ekf.nL
    length(State.Ekf.mu)
    'Warning: State vector length mismatch, ekf_hybrid_update line 85'
end
end
        