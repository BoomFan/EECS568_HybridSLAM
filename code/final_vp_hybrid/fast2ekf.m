function [ekf_mean,ekf_cov]=fast2ekf  
global pf
global com
%save=[];
ekf_mean=zeros(3+2*com,1);
ekf_cov=zeros(3+2*com,3+2*com);
%m=[];  m no longer needed
cos_theta=0;
sin_theta=0;
for j=1:length(pf)
    decide=size(pf(j).pose);
    if decide(1)==1;
        pf(j).pose=pf(j).pose';
    end
    ekf_mean(1:2,1)=pf(j).pose(1:2)/length(pf)+ekf_mean(1:2,1);
    pf(j).pose(3)=minimizedAngle(pf(j).pose(3));  %already minimized angle
    cos_theta=cos_theta+cos(pf(j).pose(3));
    sin_theta=sin_theta+sin(pf(j).pose(3));
end
ekf_mean(3,1)=atan2(1/length(pf)*sin_theta,1/length(pf)*cos_theta);  %mean of circular quantity
ekf_mean(3,1)=minimizedAngle(ekf_mean(3,1));
numb_mean=zeros(com,1);  %store the number of elements that have been averaged
numb_cov=zeros(com,com); %for covariance among landmarks purpose
numb_cov_rl=zeros(1,com); %for covariance between robot pose and landmarks

for j=1:length(pf)
    for k=1:pf(j).Mn
        n=pf(j).gama(k);  %the kth landmark in the jth particle corresponds to the nth common landmark
        ekf_mean(3+2*n-1:3+2*n,1)=ekf_mean(3+2*n-1:3+2*n,1)*numb_mean(n)+pf(j).lm(1:2,k);
        ekf_mean(3+2*n-1:3+2*n,1)=ekf_mean(3+2*n-1:3+2*n,1)/(numb_mean(n)+1);
        numb_mean(n)=numb_mean(n)+1;
    end
end    % all the common landmark mean has been calculated

%first compute covariance of robot pose:
for j=1:length(pf)
    ekf_cov(1:3,1:3)=ekf_cov(1:3,1:3)+1/length(pf)*(pf(j).pose-ekf_mean(1:3))*(pf(j).pose-ekf_mean(1:3))';
    % only the covariance among particles' robot pose are computed.
end
%Next compute covariance among landmarks
for j=1:length(pf)
    for k=1:pf(j).Mn     % for covariance, two loops of length Mn are needed
       for i=1:pf(j).Mn
m=pf(j).gama(k);
n=pf(j).gama(i);  %m corresponds to k while n corresponds to i
if m==n   %diagonal term, both P and cross covariance have to been accounted for
temp_cov=pf(j).cov(1:2,2*k-1:2*k)+(pf(j).lm(1:2,k)-ekf_mean(3+2*m-1:3+2*m,1))*(pf(j).lm(1:2,k)-ekf_mean(3+2*m-1:3+2*m,1))';
%diagonal term plus cross term
ekf_cov(3+2*m-1:3+2*m,3+2*m-1:3+2*m)=ekf_cov(3+2*m-1:3+2*m,3+2*m-1:3+2*m)*numb_cov(m,m)+temp_cov;
ekf_cov(3+2*m-1:3+2*m,3+2*m-1:3+2*m)=ekf_cov(3+2*m-1:3+2*m,3+2*m-1:3+2*m)/(numb_cov(m,m)+1);
numb_cov(m,m)=numb_cov(m,m)+1;
else   %m~=n
   temp_cov=(pf(j).lm(1:2,k)-ekf_mean(3+2*m-1:3+2*m,1))*(pf(j).lm(1:2,i)-ekf_mean(3+2*n-1:3+2*n,1))';
   ekf_cov(3+2*m-1:3+2*m,3+2*n-1:3+2*n)=ekf_cov(3+2*m-1:3+2*m,3+2*n-1:3+2*n)*numb_cov(m,n)+temp_cov;
   ekf_cov(3+2*m-1:3+2*m,3+2*n-1:3+2*n)=ekf_cov(3+2*m-1:3+2*m,3+2*n-1:3+2*n)/(numb_cov(m,n)+1);
   numb_cov(m,n)=numb_cov(m,n)+1;
end % end if
       end
    end
end

%Finally compute covariance between robot pose and landmarks 
for j=1:length(pf)
    for k=1:pf(j).Mn
        m=pf(j).gama(k);
        temp=pf(j).pose-ekf_mean(1:3,1); % for minimize angle purpose 
        temp(3)=minimizedAngle(temp(3));
        temp_cov=temp*(pf(j).lm(1:2,k)-ekf_mean(3+2*m-1:3+2*m,1))';
        ekf_cov(1:3,3+2*m-1:3+2*m)=ekf_cov(1:3,3+2*m-1:3+2*m)*numb_cov_rl(m)+temp_cov;
        ekf_cov(1:3,3+2*m-1:3+2*m)=ekf_cov(1:3,3+2*m-1:3+2*m)/(numb_cov_rl(m)+1);
        ekf_cov(3+2*m-1:3+2*m,1:3)=ekf_cov(3+2*m-1:3+2*m,1:3)*numb_cov_rl(m)+temp_cov';
        ekf_cov(3+2*m-1:3+2*m,1:3)=ekf_cov(3+2*m-1:3+2*m,1:3)/(numb_cov_rl(m)+1);
        numb_cov_rl(m)=numb_cov_rl(m)+1;
    end
end

% All the mean and covariance have been computed ,but some might be zero. Those
% should be removed: The removal criterion is if the diagonal covariance is
% zero
Ex=[];
for i=1:com
    if ekf_cov(3+2*i-1,3+2*i-1)==0 & ekf_cov(3+2*i,3+2*i)==0  %indicate a common landmark that is not in any particle
        Ex=[Ex,3+2*i-1,3+2*i];
%         % remove this landmark
%         if i~=com  %this works if this landmark is not the final one
%         ekf_mean=[ekf_mean(1:3+2*i-2,1);
%             ekf_mean(3+2*i+1:end,1)];
%         ekf_cov=[ekf_cov(1:3+2*i-2,1:3+2*i-2),  ekf_cov(1:3+2*i-2,3+2*i+1:end);
%             ekf_cov(1:3+2*i-2,3+2*i+1:end)'  ,   ekf_cov(3+2*i+1:end,3+2*i+1:end)];
%         else  %if i==com
%             ekf_mean=ekf_mean(1:3+2*com-2,1);
%             ekf_cov=ekf_cov(1:3+2*com-2,1:3+2*com-2);
%         end
    end
end
Extract=[];
for k=1:3+2*com
    if (length(find(k==Extract))==0)
        Extract=[Extract,k];
    end
end
 ekf_mean=ekf_mean(Extract);
 ekf_cov=ekf_cov(Extract,Extract);
end
    



