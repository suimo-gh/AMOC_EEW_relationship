%% create MOC data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;

modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';

expid = ''; % piControl historical ssp126 ssp245 ssp370 ssp585 u03-hos g01-hos
var = ''; % msftyz,msftmz
pathi =[''];
patho =[''];
i = []; %%% pick model %%%
% ------------------------------------------------------------------------

mons = 1:12;
rou  = 1000; %kg/m3;
for ir=[ ];
model = modellist{i};    
    clear temp
    filelist = dir(pathi); %in 'kg s-1'
    for n =1:size(filelist,1);
        temp = ncread([filelist(n).folder,'/',filelist(n).name], var);
        temp = squeeze(temp);
        temp = reshape(temp,[size(temp,1),size(temp,2),size(temp,3),12,size(temp,4)/12]);
        temp = squeeze(nanmean(temp(:,:,:,mons,:),4));
        if n ==1;
            MOC = temp;
        else
            MOC = cat(4,MOC,temp);
        end;
    end;  
MOC = MOC./rou./10^6; % convert 'kg s-1' to 'm3/s' to 'Sv'
% 3 basins 
MOC_At = squeeze(MOC(:,:,[],:)); % 1 or 2 or 3
MOC_IP = squeeze(MOC(:,:,[],:)); % 1 or 2 or 3
MOC_Gl = squeeze(MOC(:,:,[],:)); % 1 or 2 or 3
 
    file = [filelist(n).folder,'/',filelist(n).name];
    if i==22; lev = ncread(file,'olevel_bounds');lev=lev(2,:); % only for IPSL
    else lev = ncread(file,'lev'); end; % load ocean depth coordinate in m
    if i ==6 || 7; lev = lev/100; end; % for CESM models %%%%%%%

    if var == 'msftmz';     
        lat =ncread(file,'lat'); 
    elseif var == 'msftyz'; 
        lat =ncread(file,'rlat'); 
    end;
% ------------------------------------------------------------------------
%interp to standard y-z grid
levq = [0 5,10:10:190, 200:20:280, 300:100:900, 1000:200:5800];
latq = [-89.5 :1: 89.5];

clear temp; for i=1:size(MOC_At,2);for j=1:size(MOC_At,3);
    temp(:,i,j)=interp1(lat,squeeze(MOC_At(:,i,j)),latq);
end;end;MOC_At = temp;
clear temp; for i=1:size(MOC_At,1);for j=1:size(MOC_At,3);
    temp(i,:,j)=interp1(lev,squeeze(MOC_At(i,:,j)),levq);
end;end;MOC_At = temp;
save([patho,var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_At_',startyr,'_',num2str(size(MOC_At,3)),'yr_ann_interp.mat'],'MOC_At','latq','levq');

clear temp; for i=1:size(MOC_IP,2);for j=1:size(MOC_IP,3);
    temp(:,i,j)=interp1(lat,squeeze(MOC_IP(:,i,j)),latq);
end;end;MOC_IP = temp;
clear temp; for i=1:size(MOC_IP,1);for j=1:size(MOC_IP,3);
    temp(i,:,j)=interp1(lev,squeeze(MOC_IP(i,:,j)),levq);
end;end;MOC_IP = temp;
save([patho,var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_IP_',startyr,'_',num2str(size(MOC_At,3)),'yr_ann_interp.mat'],'MOC_IP','latq','levq');
% % 
clear temp; for i=1:size(MOC_Gl,2);for j=1:size(MOC_Gl,3);
    temp(:,i,j)=interp1(lat,squeeze(MOC_Gl(:,i,j)),latq);
end;end;MOC_Gl = temp;
clear temp; for i=1:size(MOC_Gl,1);for j=1:size(MOC_Gl,3);
    temp(i,:,j)=interp1(lev,squeeze(MOC_Gl(i,:,j)),levq);
end;end;MOC_Gl = temp;
save([patho,var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_Gl_',startyr,'_',num2str(size(MOC_At,3)),'yr_ann_interp.mat'],'MOC_Gl','latq','levq');

end;


%% calculate IP MOC for some models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;
modellist = {'CAS-ESM2-0','CIESM','FGOALS-f3-L','FGOALS-g3'}';

expid = ''; % piControl historical ssp126 ssp245 ssp370 ssp585 u03-hos g01-hos
var = ''; % msftyz,msftmz
pathi =[''];
i = []; %%% pick model %%%

for ir=[];
model = modellist{i}
    clear temp;
    filelist = dir(pathi); %in 'kg s-1'
load([filelist(1).folder,'/proc/',var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_Gl_',startyr,'_86yr_ann_interp.mat']); 
load([filelist(1).folder,'/proc/',var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_At_',startyr,'_86yr_ann_interp.mat']);

MOC_IP = MOC_Gl-MOC_At;
save([filelist(1).folder,'/proc/',var,'_inSv_',expid,'_run',num2str(ir),'_',model,'_MOC_IP_',startyr,'_86yr_ann_interp.mat'],'MOC_IP','latq','levq');
end;



%% Create thetao/uo/vo/wo data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;
modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';

var = ''; % thetao, uo, vo, wo
expid = ''; % piControl historical ssp126 ssp245 ssp370 ssp585 u03-hos g01-hos
pathi =[''];
patho =[''];

for im=[ ];
    model = modellist{im}
    filelist = dir(pathi); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name], var);   temp(temp>10^3) = NaN;
        if im==22;lev = ncread([filelist(1).folder,'/',filelist(1).name],'olevel'); % IPSL
        else lev = ncread([filelist(1).folder,'/',filelist(1).name], 'lev');end;
        if i ==6 || 7; lev = lev/100; end; % for CESM models %%%%%%%

        lon = ncread([filelist(1).folder,'/',filelist(1).name], 'lon');
        lat = ncread([filelist(1).folder,'/',filelist(1).name], 'lat');

% %% interp2 to standard grid
levq = [0 5,10:10:190, 200:20:280, 300:100:900, 1000:200:5800];
for i=1:size(temp,1);for j=1:size(temp,2);for y=1:size(temp,4);
    thetao(i,j,:,y)=interp1(lev,squeeze(temp(i,j,:,y)),levq);
end;end;end;

lev = levq;
save([patho,var,'_',model,'_',expid,'_run1_2080-2099_360x180x58_ann.mat'],'thetao','lev','lat','lon');

squeeze(nanmean(nanmean(thetao,1),4));
figure(im);pcolor(ans);

end;









