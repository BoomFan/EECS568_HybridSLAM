function pf=resample(pf)
% do some stuff here
pff=pf;
size_samples=size(pf,2);
sum_weight=0;
for i=1:size_samples
    sum_weight=pf(i).wt+sum_weight;
end
c(1)=pf(1).wt/sum_weight;
for i=2:size_samples
    c(i)=c(i-1)+pf(i).wt/sum_weight;
end
i=1;
u(1)=1/size_samples*rand(1);
for j=1:size_samples
    while u(j) > c(i)
        i=i+1;

    end
    pf(j)=pff(i);
    u(j+1)=u(j)+1/size_samples;
end

% for i=1:size_samples
%     pf(i).wt=1/size_samples;
% end