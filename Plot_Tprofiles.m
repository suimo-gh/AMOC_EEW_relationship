clc; clear all;

yearhist = [1980:1999];
yearfut  = [2080:2099];
pathi1 =[''];
pathi2 =[''];
expid1 = ''; 
expid2 = ''; 

modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';
% ------------------------------------------------------------------------
varname = 'thetao';
im=[];
for i=im;  
    model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    to_hist(:,:,:,i) = nanmean(thetao,4);
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    to_fut(:,:,:,i) = nanmean(thetao,4);
end;'thetao'

varname = 'uo';
im=[];
for i=im; 
    model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    uo_hist(:,:,:,i) = nanmean(thetao,4);
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    uo_fut(:,:,:,i) = nanmean(thetao,4);
end;'uo'

    latsst = ncread([filelist(1).folder,'/',filelist(1).name],'lat');
    lonsst = ncread([filelist(1).folder,'/',filelist(1).name],'lon');

clear moc_ip_hist;clear moc_ip_fut;clear moc_at_hist;clear moc_at_fut;clear moc_gl_hist;clear moc_gl_fut;
for i=1:30; i
    model = modellist{i};
        filelist = dir([pathi1,'/msft*_inSv_',expid1,'_run1_',model,'_MOC_At_ann_interp.mat']);
        load([filelist(1).folder,'/',filelist(1).name]); temp = MOC_At; temp(temp>10^10)=NaN;
    moc_at_hist(:,:,:,i) = temp(:,:,[yearhist-1850+1]);
end;clear temp;moc_at_hist=squeeze(nanmean(moc_at_hist,3));
for i=1:30; i
    model = modellist{i};
        filelist = dir([pathi2,'/msft*_inSv_',expid2,'_run1_',model,'_MOC_At_ann_interp.mat']);
        load([filelist(1).folder,'/',filelist(1).name]); temp = MOC_At; temp(temp>10^10)=NaN;
    moc_at_fut(:,:,:,i) = temp(:,:,[yearfut-2015+1]);
end;clear temp;moc_at_fut=squeeze(nanmean(moc_at_fut,3));
clear MOC_At;
 
% ------------------------------------------------------------------------
% plot clim and change - latitude mean
lonrange = [30:360-80];% tropical Indian-Pacific
% lonrange = [130:360-80]; %tropical Pacific 

hist = squeeze(nanmean(to_hist(lonrange,:,:,:),1));
fut  = squeeze(nanmean(to_fut (lonrange,:,:,:),1));

im=[];
data = nanmean(hist(:,:,im),3); mm=12; flag = 0; % clim
% data = nanmean(fut(:,:,im),3) - nanmean(hist(:,:,im),3); flag =1;mm=1; % changes

close all;figure(1); set(gcf, 'position', [0 100 400*1 160*1]); handaxes1 = axes('Position', [0.2 0.17 0.48 0.69]);hold on;set(gcf,'color','w');
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

[latp,pfullp] =meshgrid(lat,lev);fontsi = 10.5;
hold on;    [C0,h0]=contourf(latp,pfullp,data',[-50:0.1:50]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
caxis([-3 3]*mm); 

MOChist= moc_ip_hist(:,:,:); MOCfut = moc_ip_fut(:,:,:); 
data = nanmean(MOChist,3); mmm = 1;% clim
% data = nanmean(MOCfut,3) - nanmean(MOChist,3); mmm = 0.25;% changes
[latp,pfullp] =meshgrid(lat,lev(1:end));
hold on;    [C1,h1]=contour(latp,pfullp,datau',[2:2:500]*2*mmm,'linestyle','-','linecolor',rgb('bright blue'),'linewidth',0.8);  
hold on;    [C2,h2]=contour(latp,pfullp,datau',-1*[2:2:400]*2*mmm,'linestyle','--','linecolor',rgb('bright blue'),'linewidth',0.8);  
l0=plot(-90,0,'w-');

ax = gca; box on;
        ax.XTick=[-30:10:70]-0.5; ax.XLim =[-25 25]-0.5;ax.XTickLabel={'30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN'}; 
        ax.FontSize=fontsi;
         ylabel('Depth (m)','FontSize',fontsi-0.5);
        ax.YTickLabel={'0','5','10' '20','50' '100','250' '500' '1000' ,'2000' '4000'}; 
        ax.YTick=[1,5,10 20 50 100 250 500 1000 2000 4000];
        ax.YLim =[5 4000]
        set(gca,'tickdir','out');
        set(gca,'ydir','reverse','fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        ax.YScale='log';

% ------------------------------------------------------------------------
% plot clim and change - longtitude mean
latrange1 = near1(lat,-5):near1(lat,5);

hist = squeeze(mean(to_hist(:,latrange1,:,:),2)) ;
fut  = squeeze(mean(to_fut (:,latrange1,:,:),2)) ;

im=[]; 
% data = nanmean(hist(:,:,im),3); mm=12;flag =1; mm=1; % clim
data = nanmean(fut(:,:,im),3) - nanmean(hist(:,:,im),3); mm=1;flag =0; %change

close all;figure(1); set(gcf, 'position', [0 100 400*1 160*1]); handaxes1 = axes('Position', [0.2 0.25 0.5 0.71]);hold on;set(gcf,'color','w');
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

[latp,pfullp] =meshgrid(lon,lev);fontsi = 10.5;
hold on;    [C0,h0]=contourf(latp,pfullp,data',[-10:0.1:10]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
    caxis([-3 3]*1*mm); 
    c=colorbar('horizontal','position',[0.26 0.1 0.44 0.036]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-36:1:16])*mm;c.TickLength =[0 0];cbarrow

clim = nanmean(hist,3);
hold on;    [C1,h1]=contour(latp,pfullp,clim',[20 20],'linestyle','-','linecolor',addcolorplus(128),'linewidth',2);  

im=[];
latrange1 = near1(lat,-5):near1(lat,5);
uo = squeeze(nanmean(nanmean(uo_fut(:,latrange1,:,im),2),4)) - squeeze(nanmean(nanmean(uo_hist(:,latrange1,:,im),2),4)); mmm=0.2;
[latp,pfullp] =meshgrid(lon,lev(1:end));
hold on;    [C1,h1]=contour(latp,pfullp,uo',[0.05:0.05:5]*1 * mmm,'linestyle','-','linecolor',addcolorplus(1),'linewidth',0.5);  
hold on;    [C1,h1]=contour(latp,pfullp,uo',[0.05:0.05:5]*-1 * mmm,'linestyle','--','linecolor',addcolorplus(1),'linewidth',0.5);  

ax = gca; box on;
        ax.XTick=[0:30:360]; ax.XLim =[150 360-80];
        ax.XTickLabel={'0^o','30^oE','60^oE','90^oE','120^oE','150^oE','180^o','150^oW','120^oW','90^oW','60^oW','30^oW','0^o'}; 
        ax.FontSize=fontsi;
        ax.YTickLabel = {''};
        ylabel('Depth(m)','FontSize',fontsi);
        ax.YTickLabel={'0','5','10' '20','50' '100','250' '500' '1000' ,'2000' '4000'}; 
        ax.YTick=[0,5,10 20 50 100 250 500 1000 2000 4000];
        ax.YLim =[5 4000]
        set(gca,'tickdir','out');
        set(gca,'ydir','reverse','fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        ax.YScale='log';


% ------------------------------------------------------------------------
% plot regression - latitude mean
lonrange = [30:360-80];flag =1;

MOChist= moc_ip_hist; 
MOCfut = moc_ip_fut; 

im=[];
hist = squeeze(nanmean(to_hist(lonrange,:,:,im),1));
fut  = squeeze(nanmean(to_fut (lonrange,:,:,im),1)); 
amoc_max_diff=[-8.51949278573794	-8.40197415486865	-5.67319697680795	-3.80951891406736	-3.96343232444736	-15.0410973871322	-13.4270047968499	-14.9252349047024	-6.56875615789509	-5.13780915944112	-10.5574300074226	-9.27216128847374	-6.18151031875871	-4.20863883172179	-7.56060870873356	-9.10073703072974	-11.0912727105611	-8.65379637434922	-6.87719377184798	-4.59023152746667	-2.31414827840000	-4.98608210086822	-4.46318017293646	-6.79151558788571	-5.02666317164445	-5.64837601777778	-15.8516325514917	-14.1840427498667	-13.6554720280000	-7.69591656083434];
% ---------------- ---------------- ---------------- ----------------
clear reg_moc;clear tt;clear X;
for i=1:size(fut,1);for j=1:size(fut,2);
        x = (amoc_max_diff')*-1  ;  mm =50; mm = 5;%%%
        y = squeeze(fut(i,j,:)-hist(i,j,:));
        y = y([im]);
        x = x([im]);
        X=[ones(length(x),1),x];
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2) * 8; warning off;  %%% multiple by 12? 8? 10?
        tt(i,j) = stats(3);        
end;end;

close all;figure(1); set(gcf, 'position', [0 100 400*1 160*1]); handaxes1 = axes('Position', [0.2 0.17 0.48 0.69]);hold on;set(gcf,'color','w');
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);
data = reg_moc; %%%
[latp,pfullp] =meshgrid(lat,lev(1:end));fontsi = 10.5;
hold on;    [C0,h0]=contourf(latp,pfullp,data',[-30:0.01:30],'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
        caxis([-2 2]*0.2*mm); 
% c=colorbar('vertical','position',[0.93 0.28 0.015 0.6]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:2:30]*0.1*mm);cbarrow

mask=tt<0.05;stipple(latp,pfullp,mask','density',120,'color',[0 0 0]+0.5,'marker','.','markersize',5);

for i=1:size(MOCfut,1);for j=1:size(MOCfut,2);
        x = amoc_max_diff'*-1  ; %%%      
        X=[ones(length(x),1),x];
        y = squeeze(MOCfut(i,j,:)-MOChist(i,j,:));
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;
        tt(i,j) = stats(3);
end;end;
datau=reg_moc ; mmm = 0.25;
[latp,pfullp] =meshgrid(lat,lev);
hold on;    [C1,h1]=contour(latp,pfullp,datau',[2:2:500]*2*mmm,'linestyle','-','linecolor',rgb('bright blue'),'linewidth',0.8);  
hold on;    [C2,h2]=contour(latp,pfullp,datau',-1*[2:2:400]*2*mmm,'linestyle','--','linecolor',rgb('bright blue'),'linewidth',0.8);  

ax = gca; box on;
        ax.XTick=[-30:10:70]-0.5; ax.XLim =[-25 25]-0.5;ax.XTickLabel={'30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN'}; 
        ax.FontSize=fontsi;
        ylabel('Depth (m)','FontSize',fontsi-0.5);
        ax.YTickLabel={'0','5','10' '20','50' '100','250' '500' '1000' ,'2000' '4000'}; 
        ax.YTick=[0,5,10 20 50 100 250 500 1000 2000 4000];
        ax.YLim =[5 4000]
        set(gca,'tickdir','out');
        set(gca,'ydir','reverse','fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        ax.YScale='log';



% ------------------------------------------------------------------------
% plot regression - latitude mean
latrange1 = near1(lat,-5):near1(lat,5);

hist = squeeze(mean(to_hist(:,latrange1,:,:),2)) ;
fut  = squeeze(mean(to_fut (:,latrange1,:,:),2)) ;flag =1; 

amoc_max_diff=[-8.51949278573794	-8.40197415486865	-5.67319697680795	-3.80951891406736	-3.96343232444736	-15.0410973871322	-13.4270047968499	-14.9252349047024	-6.56875615789509	-5.13780915944112	-10.5574300074226	-9.27216128847374	-6.18151031875871	-4.20863883172179	-7.56060870873356	-9.10073703072974	-11.0912727105611	-8.65379637434922	-6.87719377184798	-4.59023152746667	-2.31414827840000	-4.98608210086822	-4.46318017293646	-6.79151558788571	-5.02666317164445	-5.64837601777778	-15.8516325514917	-14.1840427498667	-13.6554720280000	-7.69591656083434];
clear tt;
im=[];
for i=1:size(fut,1);for j=1:size(fut,2);
        x = (amoc_max_diff')*-1  ; %%%
        y = squeeze(fut(i,j,:)-hist(i,j,:));
        y = y([im]);
        x = x([im]);
        X=[ones(length(x),1),x]; mm = 5;
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;
        tt(i,j) = stats(3);
end;end;

close all;figure(1); set(gcf, 'position', [0 100 400*1 160*1]); handaxes1 = axes('Position', [0.2 0.25 0.5 0.71]);hold on;set(gcf,'color','w');
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);


[latp,pfullp] =meshgrid(lon,lev);fontsi = 10.5;
data = reg_moc'; %%%
hold on;    [C0,h0]=contourf(latp,pfullp,data,[-10:0.01:10]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
        caxis([-2 2]*0.2*mm); 
    c=colorbar('horizontal','position',[0.3 0.1 0.4 0.036]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-12:2:16])*0.1*mm;c.TickLength =[0 0];cbarrow

mask=tt<0.1;stipple(latp,pfullp,mask','density',180,'color',[0 0 0]+0.5,'marker','.','markersize',6);

clim = nanmean(hist,3);
hold on;    [C1,h1]=contour(latp,pfullp,clim',[20 20],'linestyle','-','linecolor',addcolorplus(128),'linewidth',2);  

% Uo-----
im=[];
hist = squeeze(mean(uo_hist(:,latrange1,:,:),2)) ;
fut  = squeeze(mean(uo_fut (:,latrange1,:,:),2)) ;
for i=1:size(fut,1);for j=1:size(fut,2);
        x = (amoc_max_diff')*-1  ; %%%
        y = squeeze(fut(i,j,:)-hist(i,j,:));
        y = y([im]);
        x = x([im]);
        X=[ones(length(x),1),x];
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;
        tt(i,j) = stats(3);
end;end;
uo = reg_moc'; %%%
hold on;    [C1,h1]=contour(latp,pfullp,uo,[0.1:0.1:5]*0.5 * mmm,'linestyle','-','linecolor',addcolorplus(1),'linewidth',0.5);  
hold on;    [C1,h1]=contour(latp,pfullp,uo,[0.1:0.1:5]*-0.5 * mmm,'linestyle','--','linecolor',addcolorplus(1),'linewidth',0.5);  

ax = gca; box on;
        ax.XTick=[0:30:360]; ax.XLim =[150 360-80];
        ax.XTickLabel={'0^o','30^oE','60^oE','90^oE','120^oE','150^oE','180^o','150^oW','120^oW','90^oW','60^oW','30^oW','0^o'}; 
        ax.FontSize=fontsi;
        ax.YTickLabel = {''};
        ylabel('Depth (m)','FontSize',fontsi);
        ax.YTickLabel={'0','5','10' '20','50' '100','250' '500' '1000' ,'2000' '4000'}; 
        ax.YTick=[0,5,10 20 50 100 250 500 1000 2000 4000];
        ax.YLim =[5 4000]
        set(gca,'tickdir','out');
        set(gca,'ydir','reverse','fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        ax.YScale='log';

