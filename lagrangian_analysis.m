%       ***************************************************
%       *  Copyright (C) 2017, Hiroshi Ashikaga, MD, PhD  *
%       *  hashika1@jhmi.edu                              *
%       *  Cardiac Arrhythmia Service                     *
%       *  Johns Hopkins University School of Medicine    *
%       *  Baltimore, Maryland, USA                       *
%       *  5/22/2017                                      *
%       ***************************************************

clear all
close all
t_end = 23810;                              % Last frame; 2.52ms/frame x 23,810 frames = 60.0012 sec

%% Information flow parameter setup

load uvo60.mat;                             % Time series of information flow
uo = uo(:,:,1:t_end);                       % Trim time series to ~60 sec
vo = vo(:,:,1:t_end);
t_int = 1.008;                              % Time window e.g. 400 frames x 2.52ms = 1.008 sec
to = 0; tmax = 0.00252*t_end;               % Frames in [sec]
ko = 0; kmax = 118.8;                       % Size in [mm] = 120 unit x 0.99 mm/unit
timespan = [to 0.00252*t6];                 % Entire time span in [sec]
domain = [ko kmax;ko kmax];                 % 2-D lattice size in [mm]
io = sqrt(uo.^2+vo.^2)/t_int;               % Magnitude of information flow  
                                            % = Information content transferred out of X per unit time [bits/sec]
avg_dist = 0.99*(4*1+4*sqrt(2))/8;          % Average distance between X and its eight Moore neighbors [mm]
vxo = avg_dist * uo./io/t_int;              % Information velocity (x component) [mm/sec] 
vyo = avg_dist * vo./io/t_int;              % Information velocity (y component) [mm/sec] 
k = linspace(ko,kmax,size(uo,1))';          % 2-D lattice size in [mm]
[k1,k2] = ndgrid(k,k);                      
time = linspace(to,tmax,size(uo,3))';       % Time in [sec]
dt = (tmax - to)/size(uo,3);                % Time step between frames in [sec]

%% Forward trajectory
xxo1 = zeros(size(vxo));                    % Forward trajectory in [mm]
xyo1 = zeros(size(vxo));                    
xxo1(:,:,1) = k1';                          % Starting points in [mm]
xyo1(:,:,1) = k2';
for frame=1:numel(time)-1                   % Integration forward in time from t0 to tmax
    
    % Liner interpolation of information velocity anywhere in 2-D lattice
    vxInterp = griddedInterpolant(k1,k2,vxo(:,:,frame)','linear');  
    vyInterp = griddedInterpolant(k1,k2,vyo(:,:,frame)','linear');
    
    % Integration of information velocity forward in time
    xxo1(:,:,frame+1) = xxo1(:,:,frame) + vxInterp(xxo1(:,:,frame)',xyo1(:,:,frame)')'*dt;
    xyo1(:,:,frame+1) = xyo1(:,:,frame) + vyInterp(xxo1(:,:,frame)',xyo1(:,:,frame)')'*dt;

end

% Show forward trajectory

figure; hold on; axis image off; set(gcf,'position',[600 350 512 512],'color',[0 0 0])
plot(xxo1(1:2:end,1:2:end,t6),xyo1(1:2:end,1:2:end,t6),'w.','markersize',40);
for m=1:2:size(xxo1,1)
    for n=1:2:size(xxo1,2)
        plot(squeeze(xxo1(m,n,1:t6)),squeeze(xyo1(m,n,1:t6)),'w','LineWidth',0.001);
    end
end
axis ij equal off; axis([11 110 11 110]);

%% Backward trajectory
xxo2 = zeros(size(vxo));                    % Backward trajectory in [mm]
xyo2 = zeros(size(vxo));
xxo2(:,:,end) = k1';                        % Ending points in [mm]
xyo2(:,:,end) = k2';
for frame=numel(time):-1:2                  % Integration backward in time from tmax to t0
    
    % Liner interpolation of information velocity anywhere in 2-D lattice
    vxInterp = griddedInterpolant(k1,k2,vxo(:,:,frame)','linear');
    vyInterp = griddedInterpolant(k1,k2,vyo(:,:,frame)','linear');
    
    % Integration of information velocity backward in time
    xxo2(:,:,frame-1) = xxo2(:,:,frame) - vxInterp(xxo2(:,:,frame)',xyo2(:,:,frame)')'*dt;
    xyo2(:,:,frame-1) = xyo2(:,:,frame) - vyInterp(xxo2(:,:,frame)',xyo2(:,:,frame)')'*dt;

end

%% Forward finite-time Lyapunov exponent (FTLE)

M = xxo1(:,:,1);                            % Reference configuration vector M,N in forward trajectory
N = xyo1(:,:,1);                            % (= 1st frame at t=0 sec)

m = xxo1(:,:,t_end);                        % Deformed configuration vector m,n in forward trajectory
n = xyo1(:,:,t_end);                        % (= Last frame at t=60 sec)

MM = [[0 M(2,:) 0];[M(:,2) M M(:,end-1)];[0 M(end-1,:) 0]]; % Padding for Neumann boundary conditions
NN = [[0 N(2,:) 0];[N(:,2) N N(:,end-1)];[0 N(end-1,:) 0]];
mm = [[0 m(2,:) 0];[m(:,2) m m(:,end-1)];[0 m(end-1,:) 0]];
nn = [[0 n(2,:) 0];[n(:,2) n n(:,end-1)];[0 n(end-1,:) 0]];

T = diff(timespan);

F11 = (mm(2:end-1,3:end)-mm(2:end-1,1:end-2))./(MM(2:end-1,3:end)-MM(2:end-1,1:end-2)); 
F12 = (mm(3:end,2:end-1)-mm(1:end-2,2:end-1))./(NN(3:end,2:end-1)-NN(1:end-2,2:end-1)); 
F21 = (nn(2:end-1,3:end)-nn(2:end-1,1:end-2))./(MM(2:end-1,3:end)-MM(2:end-1,1:end-2)); 
F22 = (nn(3:end,2:end-1)-nn(1:end-2,2:end-1))./(NN(3:end,2:end-1)-NN(1:end-2,2:end-1)); 
F11 = F11(:); F12 = F12(:); F21 = F21(:); F22 = F22(:);

for i=1:length(F11)
    F = [F11(i) F12(i);F21(i) F22(i)];      % Deformation gradient tensor F
    C = F'*F;                               % Right Cauchy-Green deformation tensor C = F'F
    C(isnan(C)==1)=0; C(isinf(C)==1)=0;     % Remove NaN and Inf
    lambda = eig(C);                        % Calculate eigenvalues of C
    lambda_m = sort(lambda(1:end));         % Sort eigenvalues in ascending order
    FTLE(i) = 1/T*log(sqrt(lambda_m(end))); % Calculate forward FTLE with largest eigenvalues
end

FTLE_forward = reshape(FTLE,[size(xxo1,1) size(xxo1,2)]); % Reshape forward FTLE

%% Backward finite-time Lyapunov exponent (FTLE)

m = xxo2(:,:,1);                            % Deformed configuration vector m,n in backward trajectory
n = xyo2(:,:,1);                            % (= Last frame at t=0 sec)

M = xxo2(:,:,t_end);                        % Reference configuration vector M,N in backward trajectory 
N = xyo2(:,:,t_end);                        % (= 1st frame at t=60 sec)

MM = [[0 M(2,:) 0];[M(:,2) M M(:,end-1)];[0 M(end-1,:) 0]]; % Padding for Neumann boundary condition
NN = [[0 N(2,:) 0];[N(:,2) N N(:,end-1)];[0 N(end-1,:) 0]];
mm = [[0 m(2,:) 0];[m(:,2) m m(:,end-1)];[0 m(end-1,:) 0]];
nn = [[0 n(2,:) 0];[n(:,2) n n(:,end-1)];[0 n(end-1,:) 0]];

T = diff(timespan);

F11 = (mm(2:end-1,3:end)-mm(2:end-1,1:end-2))./(MM(2:end-1,3:end)-MM(2:end-1,1:end-2)); 
F12 = (mm(3:end,2:end-1)-mm(1:end-2,2:end-1))./(NN(3:end,2:end-1)-NN(1:end-2,2:end-1)); 
F21 = (nn(2:end-1,3:end)-nn(2:end-1,1:end-2))./(MM(2:end-1,3:end)-MM(2:end-1,1:end-2)); 
F22 = (nn(3:end,2:end-1)-nn(1:end-2,2:end-1))./(NN(3:end,2:end-1)-NN(1:end-2,2:end-1)); 
F11 = F11(:); F12 = F12(:); F21 = F21(:); F22 = F22(:);

for i=1:length(F11)
    F = [F11(i) F12(i);F21(i) F22(i)];      % Deformation gradient tensor F
    C = F'*F;                               % Right Cauchy-Green deformation tensor C = F'F
    C(isnan(C)==1)=0; C(isinf(C)==1)=0;     % Remove NaN and Inf
    lambda = eig(C);                        % Calculate eigenvalues of C
    lambda_m = sort(lambda(1:end));         % Sort eigenvalues in ascending order
    FTLE(i) = 1/T*log(sqrt(lambda_m(end))); % Calculate backward FTLE with largest eigenvalues
end

FTLE_backward = reshape(FTLE,[size(xxo1,1) size(xxo1,2)]); % Reshape backward FTLE

save FTLE_forward FTLE_backward

%% Lagrangian coherent structures (LCS)

figure; 
ax1 = axes; h1 = surf(FTLE_forward,'edgealpha',0); view(2); alpha(h1,'z');
ax2 = axes; h2 = surf(FTLE_backward,'edgealpha',0); view(2); alpha(h2,'z');
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'mporange');                   % Orange = repelling LCS 
colormap(ax2,'mppink');                     % Pink = attracting LCS
caxis(ax1,[0 0.04]);
caxis(ax2,[0 0.03]);
axis(ax1,'tight','equal','off','ij'); axis(ax1,[11 110 11 110])
axis(ax2,'tight','equal','off','ij'); axis(ax2,[11 110 11 110])
set(gcf,'position',[600 350 512 512],'color',[0 0 0])
hold on; 
a1 = k1'; a2 = k2';
quiver(a1(1:4:end,1:4:end),a2(1:4:end,1:4:end),...
    xxo1(1:4:end,1:4:end,t6)-xxo1(1:4:end,1:4:end,1),...
    xyo1(1:4:end,1:4:end,t6)-xyo1(1:4:end,1:4:end,1),2,'color',[1 1 1]);
axis('tight','equal','off','ij'); axis([11 110 11 110])
axis('tight','equal','off','ij'); axis([11 110 11 110])