
clc; clear all;

yearhist = [1980:1999];
yearfut  = [2080:2099];
pathi1 =[''];
pathi2 =[''];
expid1 = ''; 
expid2 = ''; 

modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';
% ------------------------------------------------------------------------
varname = 'tos';
for i=1:30; i
    model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
    ssthist(:,:,:,i) = temp(:,:,:);
end;clear temp;ssthist=squeeze(nanmean(ssthist,3));
for i=1:30; i
    model = modellist{i};
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
    sstfut(:,:,:,i) = temp(:,:,:);
end;clear temp;sstfut=squeeze(nanmean(sstfut,3));
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

%%load colormap
clear co;
co(:,1) = addcolorplus(268);  co(:,2) = addcolorplus(270); co(:,3) = addcolorplus(264); co(:,4) = addcolorplus(266);  co(:,5) = addcolorplus(260);  co(:,6) = addcolorplus(262); 
co(:,7) = addcolorplus(256);  co(:,8) = addcolorplus(258); co(:,9) = addcolorplus(252); co(:,10) = addcolorplus(254); co(:,11) = addcolorplus(249); co(:,12) = addcolorplus(251); 

% ------------------------------------------------------------------------
clear amoc;
    latrange = near1(latq,0):near1(latq,65); 
    levrange = near1(levq,500):near1(levq,2000); 
    for i=1:size(moc_at_fut,3); 
        amoc_max_hist(i) = max(max(moc_at_hist(latrange,levrange,i)));
        amoc_max_fut(i)  = max(max(moc_at_fut (latrange,levrange,i)));
    end;
        amoc_max_diff = amoc_max_fut-amoc_max_hist;
% ------------------------------------------------------------------------
lonrange = near1(lonsst,170):near1(lonsst,360-80); latrange = near1(latsst,-5):near1(latsst,5); 
eqP   = squeeze(nanmean(nanmean(sstfut(lonrange,latrange,:),1),2)) - squeeze(nanmean(nanmean(ssthist(lonrange,latrange,:),1),2));

lonrange = near1(lonsst,130):near1(lonsst,360-80); latrange = [near1(latsst,-20):near1(latsst,20)]; % tropical P
% lonrange = near1(lonsst,30):near1(lonsst,360-80); latrange = [near1(latsst,-20):near1(latsst,20)];% tropical IP
% lonrange = near1(lonsst,0):near1(lonsst,360); latrange = [near1(latsst,-20):near1(latsst,20)]; % tropics
tropicP   = squeeze(nanmean(nanmean(sstfut (lonrange,latrange,:),1),2))-squeeze(nanmean(nanmean(ssthist (lonrange,latrange,:),1),2));
eew = eqP-tropicP;
% ------------------------------------------------------------------------
close all;figure(1); set(gcf, 'position', [0 100 360*1 360*1]); handaxes1 = axes('Position', [0.2 0.36 0.5 0.5]);hold on;set(gcf,'color','w');fontsi=11;
x = amoc_max_diff; 
y = eew'; 
for i=1:12;
    plot(x(i),y(i),'Marker','s','markersi',10,'markerfacecolor','none','markeredgecolor',co(:,i),'LineStyle','none','linewi',2);
end;
for i=13:24;'none'
    plot(x(i),y(i),'Marker','^','markersi',10,'markerfacecolor','none','markeredgecolor',co(:,i-12),'LineStyle','none','linewi',2);
end;
for i=25:30;
    plot(x(i),y(i),'Marker','o','markersi',10,'markerfacecolor','none','markeredgecolor',co(:,i-24),'LineStyle','none','linewi',2);
end;

X=[ones(size(x));x];
[b,bint,r,rint,stats] = regress(y',X') ;
plot(x,x*b(2)+b(1),'color','k','LineStyle','-','linewi',1);

ax = gca; box on;
        ax.FontSize=fontsi;
        xlabel('\Delta AMOC (Sv)','FontSize',fontsi);ax.XLim =[-19,-0];ax.XTick=[-24:6:0];
        ylabel('El Niño-like warming (°C)','FontSize',fontsi);ax.YLim =[0,0.8];ax.YTick=[0.1:0.2:1];
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        set(gca,'tickdir','out');

% ------------------------------------------------------------------------
[a,b] = corr(x',y') 
for i=1:20000;
    loc = randi([1, 30], 30-5, 1);
    [aa(i),bb(i)] = corr(x(loc)',y(loc)') ;
end;
[prctile(aa,90),prctile(aa,10)]


