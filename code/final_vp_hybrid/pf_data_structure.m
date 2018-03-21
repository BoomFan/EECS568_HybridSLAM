function pf=pf_data_structure(N,x,y,theta)
  % N is the # of particles; x,y,theta are initial pose
for p=1:N
    pf(p).wt=1/N;   %initial weight
    pf(p).pose=[x;y;theta];   %initial belief of robot pose 3*1 vecotr
    pf(p).lm=[];   %landmark x,y coordinate, 2*Mn matrix, where Mn is the # of landmarks 
    %in the p_th particle, lm augments with new observation.
    pf(p).cov=[];  %%landmark covariance, 2*2Mn matrix, where Mn is the # of landmarks 
    %in the p_th particle, cov augments with new observation.
    pf(p).Mn=0;   % # of landmarks in the nth particle
    pf(p).gama=[]; % 1*Mn vector, 'pf(p).gama(i)=n' indicates i_th landmark in p_th particle
    %corresponds to n_th landmark in the environment
    pf(p).i=[];  %1*Mn vector, counter for each landmark,see p463 FastSLAM2.0 line 25,32,39,41 and 42
    pf(p).temp=[];
    pf(p).c_bar=[];
end 
end