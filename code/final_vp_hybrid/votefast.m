function I=votefast(n,w,spur)
% this is the function to set an consensus agreement to every particle
global pf
global com;
%
%
for i=1:length(pf)
    if spur(i)==0
    A(i)=pf(i).temp(n);
    end
end
total_weight=zeros(com+1,1);
if exist('A')
for j=1:com+1
    h=find(A==j);
    for m=1:length(h)
    total_weight(j)=total_weight(j)+w(h(m));
    end
end
end
%
[Y,I] = max(total_weight);
if I==com+1
    com=com+1;
end
%%%%%
