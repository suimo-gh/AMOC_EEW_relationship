clc; clear all;

yearhist = [1980:1999];
yearfut  = [2080:2099];

levrange = near1(lev,2000):near1(lev,6500);

pathi1 =[''];
pathi2 =[''];
expid1 = ''; 
expid2 = ''; 
modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';
% ------------------------------------------------------------------------
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

varname = 'vo';
im=[];
for i=im; 
    model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    vo_hist(:,:,:,i) = nanmean(thetao,4);
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    vo_fut(:,:,:,i) = nanmean(thetao,4);
end;'vo'

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
 

r = 6400000; pi = 3.14;
dlev = diff(lev,1,2); dlev(end+1)=dlev(end);
dx = 2*pi*r/360 * grid;
A = dlev.*dx;
for i=1:size(uo_fut,1), for j=1:size(uo_fut,2), for m=1:30;
    uo_fut (i,j,:,m) = squeeze(uo_fut (i,j,:,m)).*A' ./10^6;
    uo_hist(i,j,:,m) = squeeze(uo_hist(i,j,:,m)).*A' ./10^6;
end;end;end; % u(m/s) * A(m2) / 10^6 = Sv   --- CMIP6

%%load colormap
clear co;
co(:,1) = addcolorplus(268);  co(:,2) = addcolorplus(270); co(:,3) = addcolorplus(264); co(:,4) = addcolorplus(266);  co(:,5) = addcolorplus(260);  co(:,6) = addcolorplus(262); 
co(:,7) = addcolorplus(256);  co(:,8) = addcolorplus(258); co(:,9) = addcolorplus(252); co(:,10) = addcolorplus(254); co(:,11) = addcolorplus(249); co(:,12) = addcolorplus(251); 

% ------------------------------------------------------------------------
% plot change pattern for zonal transport
hist = squeeze(nansum(uo_hist(:,:,levrange,:),3));
fut  = squeeze(nansum(uo_fut (:,:,levrange,:),3)); 
reg_moc = nanmean(fut,3) - nanmean(hist,3); mm=0.8;%%%

close all;figure(1); set(gcf, 'position', [0 100 450*1 240*1]); handaxes1 = axes('Position', [0.1 0.15 0.8 0.8]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[-90 360-90],'lat',[-80 80]); fontsi = 10;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),ones(2,3),colorfinal(16:end,:));colormap(colorfinal);

data = reg_moc'; 
data = cat(2,data,data);
[latp,pfullp] =meshgrid(cat(1,lon-360,lon),lat);
hold on;    [C0,h0]=m_contourf(latp,pfullp,data,[-5:0.05:5]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
   caxis([-2 2]*0.5*mm); 
c=colorbar('vertical','position',[0.92 0.2 0.016 0.6]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:1:30]*0.5*mm);cbarrow

lonx = cat(1,lon-360,lon);
kuang=data ;kuang(:)=0;
    lonrange = near1(lonx,20):near1(lonx,20); 
    latrange = near1(lat,-50):near1(lat,-30); 
kuang(latrange,lonrange)=100;
hold on;    [C1,h1]=m_contour(latp,pfullp,kuang,[100 100],'linestyle','-','linecolor',rgb('black'),'linewidth',2.5);  

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.8);  
m_coast('linewidth',1,'color',[0 0 0]+0.8);  
m_grid('box','on','linewidth',1,'linest','none','xtick',([-90:60:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

data = fut - hist;
data = cat(1,data,data);
SO_x_sv= squeeze(nansum(nansum(data(lonrange,latrange,:),1),2));

% ------------------------------------------------------------------------
% plot SO vs AMOC
amoc_max_diff=[-8.51949278573794	-8.40197415486865	-5.67319697680795	-3.80951891406736	-3.96343232444736	-15.0410973871322	-13.4270047968499	-14.9252349047024	-6.56875615789509	-5.13780915944112	-10.5574300074226	-9.27216128847374	-6.18151031875871	-4.20863883172179	-7.56060870873356	-9.10073703072974	-11.0912727105611	-8.65379637434922	-6.87719377184798	-4.59023152746667	-2.31414827840000	-4.98608210086822	-4.46318017293646	-6.79151558788571	-5.02666317164445	-5.64837601777778	-15.8516325514917	-14.1840427498667	-13.6554720280000	-7.69591656083434];
x = amoc_max_diff; 
y = So_x_sv'; 

im = [];
loc = setdiff(1:30,im);
y(loc)=NaN; 

close all;figure(1); set(gcf, 'position', [0 100 360*1 360*1]); handaxes1 = axes('Position', [0.2 0.36 0.5 0.5]);hold on;set(gcf,'color','w');fontsi=10;
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
        xlabel('\Delta AMOC (Sv)','FontSize',fontsi);ax.XLim =[-20,-0];
        ylabel('\Delta Zonal transport (Sv)','FontSize',fontsi);ax.YLim =[-10,0];
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        et(gca,'tickdir','out');

[a,b] = corr(x(loc)',y(loc)') 
for i=1:15000;
    loc = randi([1, 30], 25, 1);
    [aa(i),bb(i)] = corr(x(loc)',y(loc)') ;
end;
[prctile(aa,90),prctile(aa,10)]

% ------------------------------------------------------------------------
% plot change pattern for zonal transport
for i=1:size(vo_fut,1), for j=1:size(vo_fut,2), for m=1:30;
    vo_fut (i,j,:,m) = squeeze(vo_fut (i,j,:,m)).*A' ./10^6;
    vo_hist(i,j,:,m) = squeeze(vo_hist(i,j,:,m)).*A' ./10^6;
end;end;end;% u(cm/s) * A(m2) /100 / 10^6 = Sv --- CESM1

hist = squeeze(nansum(vo_hist(:,:,levrange,:),3));
fut  = squeeze(nansum(vo_fut (:,:,levrange,:),3)); 

im = [];
reg_moc = nanmean(fut(:,:,im),3) - nanmean(hist(:,:,im),3); mm = 0.6;

close all;figure(1); set(gcf, 'position', [0 100 450*1 240*1]); 
handaxes1 = axes('Position', [0.1 0.15 0.8 0.8]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[-90 360-90],'lat',[-80 80]); 
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),ones(2,3),colorfinal(16:end,:));colormap(colorfinal);

data = reg_moc';
data = cat(2,data,data);

[latp,pfullp] =meshgrid(cat(1,lon-360,lon),lat);fontsi = 10;
hold on;    [C0,h0]=m_contourf(latp,pfullp,data,[-50:0.05:50]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.5);  
        caxis([-2 2]*0.5*mm); 

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.8);  
m_coast('linewidth',1,'color',[0 0 0]+0.8);  

data = fut - hist;
data = cat(1,data,data);

lonx = cat(1,lon-360,lon);kuang=data ;kuang(:)=0;
    lonrange = near1(lonx,40):near1(lonx,70); 
    latrange = near1(lat,-20):near1(lat,-20); 
    I_y_sv= squeeze(nansum(nansum(data(lonrange,latrange,:),1),2))'; 
kuang(latrange,lonrange)=100;
hold on;    [C1,h1]=m_contour(latp,pfullp,kuang,[100 100],'linestyle','-','linecolor',rgb('black'),'linewidth',2);  
    lonrange = near1(lonx,170):near1(lonx,180+20); 
    latrange = near1(lat,-20):near1(lat,-20); 
    P_y_sv= squeeze(nansum(nansum(data(lonrange,latrange,:),1),2))'; 
kuang(latrange,lonrange)=100;
hold on;    [C1,h1]=m_contour(latp,pfullp,kuang,[100 100],'linestyle','-','linecolor',rgb('black'),'linewidth',2);  
    lonrange = near1(lonx,-40):near1(lonx,-10); 
    latrange = near1(lat,-20):near1(lat,-20); 
    At_y_sv= squeeze(nansum(nansum(data(lonrange,latrange,:),1),2))'; 
kuang(latrange,lonrange)=100;
hold on;    [C1,h1]=m_contour(latp,pfullp,kuang,[100 100],'linestyle','-','linecolor',rgb('black'),'linewidth',2);  

m_grid('box','on','linewidth',1,'linest','none','xtick',([-90:60:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])


% ------------------------------------------------------------------------
% plot IP vs AMOC
x = amoc_max_diff;
y = P_y_sv + I_y_sv;

im = [];
loc = setdiff(1:30,im);
y(loc)=NaN; 

close all;figure(1); set(gcf, 'position', [0 100 360*1 360*1]); handaxes1 = axes('Position', [0.2 0.36 0.5 0.5]);hold on;set(gcf,'color','w');fontsi=10;
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
        xlabel('\Delta AMOC (Sv)','FontSize',fontsi);ax.XLim =[-20,-0];
        ylabel('\Delta Meridional transport (Sv)','FontSize',fontsi);ax.YLim =[-10,0];
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
                set(gca,'tickdir','out');

[a,b] = corr(x(loc)',y(loc)') 
for i=1:20000;
    loc = randi([1, 30], 25, 1);
    [aa(i),bb(i)] = corr(x(loc)',y(loc)') ;
end;
[prctile(aa,90),prctile(aa,10)]

% ------------------------------------------------------------------------
% plot At vs AMOC
x = amoc_max_diff;
y = At_y_sv;

im = [];
loc = setdiff(1:30,im);
y(loc)=NaN; 

close all;figure(1); set(gcf, 'position', [0 100 360*1 360*1]); handaxes1 = axes('Position', [0.2 0.36 0.5 0.5]);hold on;set(gcf,'color','w');fontsi=10;
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
        xlabel('\Delta AMOC (Sv)','FontSize',fontsi);ax.XLim =[-18,-0];ax.XTick=[-18:6:0];
        ylabel('\Delta Meridional transport (Sv)','FontSize',fontsi);ax.YLim =[0,10];
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.6;
        set(gca,'tickdir','out');

[a,b] = corr(x(loc)',y(loc)') 
for i=1:20000;
    loc = randi([1, 30], 25, 1);
    [aa(i),bb(i)] = corr(x(loc)',y(loc)') ;
end;
[prctile(aa,90),prctile(aa,10)]
