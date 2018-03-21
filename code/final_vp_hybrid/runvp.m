function runvp(nSteps,pauseLen, makeVideo)

global Param;
global State;
global Data;
global com;
global pf;
com=0;

t1=clock
if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

Data = load_vp_si();

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================

State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);

if ~exist('makeVideo','var') || isempty(makeVideo)
    makeVideo = false;
end
if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end


global AAr;
AAr = [0:360]*pi/360;

Np=50;  %#of particles
% pf=pf_data_structure(Np,0, 0, 0);
pf=pf_data_structure(Np,Data.Gps.x(2), Data.Gps.y(2), 36*pi/180);%for fastslam use
ini_pi=[1];

figure(1); clf;
axis equal;
a=0;
ci = ones(1,Np); % control index
cc=1; %Gps index
kk=2;
vv=600;
t = min(Data.Laser.time(1), Data.Control.time(1));
predicttime=zeros(1,min(nSteps, length(Data.Laser.time)));
for k=1:min(nSteps, length(Data.Laser.time))

    k
    tic;

    for m=1:length(pf)
        Rt=0;
        while (Data.Control.time(ci(m)) < Data.Laser.time(k))
            % control available
            dt = Data.Control.time(ci(m)) - t;
            t = Data.Control.time(ci(m));
            u = [Data.Control.ve(ci(m)), Data.Control.alpha(ci(m))]';
            [ pf(m),Rt ] = predict_vp1( pf(m),u,dt,Rt);
            ci(m) = ci(m)+1;
        end
        
        % observation available
        dt = Data.Laser.time(k) - t;
        t = Data.Laser.time(k);
        z = detectTreesI16(Data.Laser.ranges(k,:));
        zz=z(1:2,:);
        zz(2,:)=minimizedAngle(z(2,:)-pi/2);
        u = [Data.Control.ve(ci(m)), Data.Control.alpha(ci(m))]';
        pftemp=pf(m);
        [ pf(m),Rt ] = predict_vp1( pf(m),u,dt,Rt);
        pf(m).wt=1;
        
%         [ new,pf(m),chat,p0,zt_zhat,Rt ] = predictvp( pf(m),u, dt,zz(:,1),Rt);       
%         if dt==0 % To prevent white noise in robot position when no action is taken
%             pf(m)=pftemp';
%         end
% %         if spur==0
%             [ w(m),pf(m) ] = update_vp( new,pf(m),chat,p0,zz(:,1),zt_zhat,Rt,1);
%             pf(m).c_bar(1)=chat;
%             pf(m).wt=w(m);
% %         end
    end
    
%     I(1)=votefast(1,w);
%     for i=1:length(pf)
%         pf(i).gama(pf(i).c_bar(1))=I(1);
%         for j=1:pf(i).Mn
%             A(j)=pf(i).gama(j);
%         end
%         GG=find(A(:)==pf(i).temp(1));
%         for f=1:length(GG)
%             pf(i).gama(GG(f))=I(1);
%         end
%         A=[];
%     end
%     
    p0=0.03;
    for n=1:size(z,2)
        for m=1:length(pf)
            [chat,new,pf(m),zt_zhat,spur(m)]=newda_nn(zz(:,n),pf(m),Rt);
            if spur(m)==0
                [w(m),pf(m)] = update_vp( new,pf(m),chat,p0,zz(:,n),zt_zhat,Rt,n )  ;
                pf(m).c_bar(n)=chat;
                pf(m).wt=pf(m).wt*w(m);
            end
        end
        
        
        I(n)=votefast(n,w,spur);
        for i=1:length(pf)
            if spur(i)==0
                pf(i).gama(pf(i).c_bar(n))=I(n);
                for j=1:pf(i).Mn
                    A(j)=pf(i).gama(j);
                end
                GG=find(A(:)==pf(i).temp(n));
                for f=1:length(GG)
                    pf(i).gama(GG(f))=I(n);
                end
                A=[];
            end
        end
    end
    
    
    while (Data.Gps.time(cc)<t)
        cc=cc+1;
    end
    while (Data.Gps.time(cc)<t)
        cc=cc+1;
    end
    pf2=pf;
    pf=resample(pf);
    
    
    for i=1:length(pf)
        pf(i).temp=[];
        pf(i).weight=[];
        pf(i).c_bar=[];
    end
    %
    
    if mod(k,vv)==0
        figure(kk);
        [ekf_mean,ekf_cov]=fast2ekf;
       
        
        % transform to truely gloabal landmark
        %          POS_mean(1:2,kk)=ekf_mean(1:2,1);
        %          POS_mean(6,kk)=ekf_mean(3,1);
        %
        % initialize the matrix
        cov_l_m=zeros(length(ekf_mean)-3,length(ekf_mean)-3);
        local_lm_mean=zeros(2,(length(ekf_mean)-3)/2);
        %prepare for head2tail
         
        X=zeros(6*(length(ekf_mean)-3)/2,1);
        J=zeros(3*(length(ekf_mean)-3)/2,6);
        
        
        [xx,cov_l_m]=head2tail_local2global(ekf_mean,ekf_cov);

        %         cov_l_m=cov_l_m*3;
         plotcov2d(xx(3+2*State.Ekf.nL+1),xx(3+2*State.Ekf.nL+2),cov_l_m(3+2*State.Ekf.nL+1:3+2*State.Ekf.nL+2,3+2*State.Ekf.nL+1:3+2*State.Ekf.nL+2), 'b', 0, 'b', 1, 3);
        hold on
        for kf=1:(length(xx)-6)/2-State.Ekf.nL
            plotcov2d(xx(6+2*State.Ekf.nL+2*kf-1), xx(6+2*State.Ekf.nL+2*kf), cov_l_m(6+2*State.Ekf.nL+2*kf-1:6+2*State.Ekf.nL+2*kf,6+2*State.Ekf.nL+2*kf-1:6+2*State.Ekf.nL+2*kf), 'r', 0, 'b', 1, 3);
        end
        axis equal
        Li=da_dgnn_modified_for_hybrid(xx(7+2*State.Ekf.nL:end), cov_l_m(7+2*State.Ekf.nL:end,7+2*State.Ekf.nL:end));%is it necessary to consider the cross correlation between global and local landmarks?
        Ekf_Hybrid_update(xx(4+2*State.Ekf.nL:end),cov_l_m,Li);
        

        %        State.Ekf.Sigma(1:3,1:3)=State.Ekf.Sigma(1:3,1:3)*1.5;
        %         max(abs(xx))
        %         max(abs(State.Ekf.mu))
        
        %
        %add code here
        figure(100+kk);
        plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'b', 0, 'b', 1, 3);
        hold on
        %         for i=1:State.Ekf.nL
        %             plotcov2d( xx(3+2*i-1), xx(3+2*i), cov_l_m(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 'r', 0, 'b', 1, 3);
        %         end
        for i=1:State.Ekf.nL
            plotcov2d( State.Ekf.mu(3+2*i-1), State.Ekf.mu(3+2*i), State.Ekf.Sigma(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 'r', 0, 'b', 1, 3);
        end
        hold off
        axis equal
        %         for i=1:(length(ekf_mean)-3)/2
        %             plotcov2d( xx(6+2*State.Ekf.nL+2*i-1), xx(6+2*State.Ekf.nL+2*i), cov_l_m(6+2*State.Ekf.nL+2*i-1:6+2*State.Ekf.nL+2*i,6+2*State.Ekf.nL+2*i-1:6+2*State.Ekf.nL+2*i), 'r', 0, 'b', 1, 3);
        %         end
        %
        kk=kk+1;
        %get new particles here
        pf=pf_data_structure(Np,0,0,0);
        com=0;
    end
    steptime(k)=toc;    
%     steptime(k)=steptime(k)/5;
    
    %----- calculate an estimate of current pose begin ------
    est_pose_x=0;
    est_pose_y=0;
    est_pose_theta_s=0;
    est_pose_theta_c=0;
    for m=1:Np
        est_pose_x=pf(m).pose(1)+est_pose_x;
        est_pose_y=pf(m).pose(2)+est_pose_y;
        est_pose_theta_s=sin(pf(m).pose(3))+est_pose_theta_s;
        est_pose_theta_c=cos(pf(m).pose(3))+est_pose_theta_c;
    end
    est_pose_x=est_pose_x/Np;
    est_pose_y=est_pose_y/Np;
%     est_pose_x=State.Ekf.mu(1);
%     est_pose_y=State.Ekf.mu(2);
    est_pose_theta_s=est_pose_theta_s/Np;
    est_pose_theta_c=est_pose_theta_c/Np;
    est_pose_theta=atan2(est_pose_theta_s,est_pose_theta_c);
    trajectory(1,k)=est_pose_x;
    trajectory(2,k)=est_pose_y;
    true_traj(1,k)=Data.Gps.x(cc);
    true_traj(2,k)=Data.Gps.y(cc);
%     %----- calculate an estimate of current pose end ------
%         doGraphics(z,cc,pf,Np,est_pose_x,est_pose_y,est_pose_theta);
%         axis([-80 0 -80 10])
%        drawnow;
       
       if k==16 || k==50 ||k==100 ||k==300 || k==500 || k==800 || k==1000 || k==1300 || k==1500 || k==1800 || k==2000
           plotfinal(State,trajectory,true_traj,steptime,k,Np)
           now_time=etime(clock,t1)
       end       
       
       
       
       
       
       
       
       
       
    if  pauseLen > 0
        pause(pauseLen);
    end
    if makeVideo
        F = getframe(gcf);
        switch votype
            case 'avifile'
                vo = addframe(vo, F);
            case 'VideoWriter'
                writeVideo(vo, F);
            otherwise
                error('unrecognized votype');
        end
    end
        
end
if makeVideo
    fprintf('Writing video...');
    switch votype
      case 'avifile'
        vo = close(vo);
      case 'VideoWriter'
        close(vo);
      otherwise
        error('unrecognized votype');
    end
    fprintf('done\n');
end

total_time=etime(clock,t1);
end

function plotfinal(State,trajectory,true_traj,steptime,k,Np)
%%%% plot the final map %%%%
figure
hold on
plot(trajectory(1,:),trajectory(2,:),'b-');
plot(true_traj(1,:),true_traj(2,:),'g-');
% for i=1:Np
%     plot(pf(i).pose(1),pf(i).pose(2),'.')
%     hold on
%     for j=1:pf(i).Mn
%         plotcov2d( pf(i).lm(1,j), pf(i).lm(2,j), pf(i).cov(1:2,2*j-1:2*j), 'b', 0, 'b', 1, 3);
%     end
% end
for i=1:State.Ekf.nL
    plotcov2d( State.Ekf.mu(2*i+1), State.Ekf.mu(2*i+2), State.Ekf.Sigma(2*i+1:2*i+2,2*i+1:2*i+2), 'b', 0, 'b', 1, 3);
end
grid on
axis equal
legend('estimate trajectory','GPS trajectory');
h_title=title({['Final map'],...
    [ 'Step=',num2str(k),', Num of particle= ',num2str(Np)]});
set(h_title,'FontSize',18);
    print(gcf,'-djpeg' ,['hybridslam_finalmap_step',num2str(k),'_Np',num2str(Np),'_p0_003_vv600..jpeg'],'-r400')
figure
hold on
p1=plot(1:k,steptime(1:k),'b-','linewidth',2);
grid on
h_title=title({['plot of CPU Time at each step'],...
    [ 'Step=',num2str(k),', Num of particle= ',num2str(Np)]});
set(h_title,'FontSize',16);
legend([p1],'time for each step','location','best')
    print(gcf,'-djpeg' ,['hybridslam_time_step',num2str(k),'_Np',num2str(Np),'_p0_003_vv600.jpeg'],'-r400')
end

%==========================================================================
function doGraphics(z,cc,pf,Np,est_pose_x,est_pose_y,est_pose_theta)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!
global Data;
figure(1)
%---------------------plot of robot---------------
%plotcov2d( pf.(1).pose(1), pf.(1).pose(2), State.Ekf.Sigma(1:2,1:2), 'r', 0, 'b', 1, 3);
for i=1:Np
    plot(pf(i).pose(1),pf(i).pose(2),'.')
    hold on
    for k=1:pf(i).Mn
        plotcov2d( pf(i).lm(1,k), pf(i).lm(2,k), pf(i).cov(1:2,2*k-1:2*k), 'b', 0, 'b', 1, 3);
    end
end


% restrict view to a bounding box around the current pose
BB=20;
axis([[-BB,BB]+est_pose_x, [-BB,BB]+est_pose_y]);


% project raw sensor detections in global frame using estimate pose
xr = est_pose_x;
yr = est_pose_y;
tr = est_pose_theta;
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end



% plot Gps 
plotbot(Data.Gps.x(cc),Data.Gps.y(cc), est_pose_theta, 'black', 1, 'red', 1);


hold off;
end