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

%% Calculate Shannon entropy

load bi60.mat;                              % Binary time series 
bmt = bmt(:,:,1:t_end);                     % Trim time series to ~60 sec
[rr,cc,ff] = size(bmt);
bmt = reshape(bmt,[rr*cc ff]);              % Reshape time series matrix for computation
b = 2;                                      % Bin number for histogram (binary: 0 or 1)
for elem=1:size(bmt,1)
    fprintf('Working on element No. %d ...\n',elem);
    Sm = bmt(elem,:); 
    hm = histcounts(Sm,b);                  
    pSm = hm/sum(hm);                       % Probability distribution 
    Hm(elem,:) = sum(-(pSm(pSm>0).*(log2(pSm(pSm>0))))); % Shannon entropy (in bits)
end
He = reshape(Hm,[rr cc]);                   % Reshape time series matrix back to original shape

% Show Shannon entropy

figure; hold on; axis image off; set(gcf,'position',[600 350 512 512],'color',[1 1 1])
imagesc(He); axis ij equal off; axis([11 110 11 110]); colormap(jet); C = He(:); 
caxis([median(C)-3*std(C) median(C)+3*std(C)]);

clear bmt

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

%% Instantaneous information flow

a1 = k1'; a2 = k2';
figure; hold on; axis image off; set(gcf,'position',[600 350 512 512],'color',[1 1 1])

% Streamlines at 1st time frame (e.g. time = 0 sec)
h1 = streamslice(a1,a2,vxo(:,:,1),vyo(:,:,1),0.5,'Method','cubic'); set(h1,'LineWidth',2,'Color','r','LineStyle','-');

% Streamlines at last time frame (e.g. time = 60 sec)
h2 = streamslice(a1,a2,vxo(:,:,t6),vyo(:,:,t6),0.5,'Method','cubic'); set(h2,'LineWidth',2,'Color','k','LineStyle','-')

axis ij equal off; axis([11 110 11 110]);

%% Total information flow over time

Svx = sum(uo(:,:,1:t_end)*dt/t_int,3);          % Vector sum of instantaneous information flow over time in [bits]
Svy = sum(vo(:,:,1:t_end)*dt/t_int,3);
Svm = sqrt(Svx.^2+Svy.^2);                      % Magnitude sum of instantaneous information flow over time in [bits]           
a1 = k1'; a2 = k2';
figure; hold on; axis image off; set(gcf,'position',[600 350 512 512],'color',[1 1 1]);
imagesc(Svm); quiver(a1(1:4:end,1:4:end),a2(1:4:end,1:4:end),Svx(1:4:end,1:4:end),Svy(1:4:end,1:4:end),2,'k');
axis ij equal off; axis([11 110 11 110]); colormap(jet); C = Svm(:); 
caxis([median(C)-3*std(C) median(C)+3*std(C)]); 