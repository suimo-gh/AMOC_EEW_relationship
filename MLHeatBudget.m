clc; clear all;

yearhist = [1980:1999];
yearfut  = [2080:2099];
pathi1 =[''];
pathi2 =[''];
expid1 = ''; 
expid2 = ''; 

modellist = {'ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2','CESM2-WACCM','CIESM','CMCC-ESM2','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-CC','FGOALS-f3-L','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-2-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'}';
% ------------------------------------------------------------------------
varname = 'hfss';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
hfss_hist = hist; 
hfss_fut  = fut; 

varname = 'hfls';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
hfls_hist = hist; 
hfls_fut  = fut; 

varname = 'rlds';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
rlds_hist = hist; 
rlds_fut  = fut; 

varname = 'rlus';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
rlus_hist = hist; 
rlus_fut  = fut; 

varname = 'rsus';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
rsus_hist = hist; 
rsus_fut  = fut; 

varname = 'rsds';
for i=[1:30]; i
model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        hist(:,:,:,i) = temp;
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],varname);
        fut(:,:,:,i) = temp;
end;
rsds_hist = hist; 
rsds_fut  = fut; 
clear hist;clear fut; clear temp
 
    latsst = ncread([filelist(1).folder,'/',filelist(1).name],'lat');
    lonsst = ncread([filelist(1).folder,'/',filelist(1).name],'lon');

% ------------------------------------------------------------------------
% plot change patterns for each terms 
close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.36 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 9.5;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

tauuhist = hfss_hist * -1;tauufut = hfss_fut *-1; %% sensitive
% tauuhist = hfls_hist * -1;tauufut = hfls_fut *-1; %% latent
% tauuhist = rlds_hist - rlus_hist;tauufut = rlds_fut - rlus_fut; % longwave
% tauuhist = rsds_hist - rsus_hist;tauufut = rsds_fut - rsus_fut; % shortwave
% tauuhist= (hfss_hist+hfls_hist)*-1 +(rlds_hist-rlus_hist)+(rsds_hist-rsus_hist);tauuhist= 1*tauuhist;%Qnet
% tauufut = (hfss_fut+hfls_fut) * -1 +(rlds_fut - rlus_fut)+(rsds_fut - rsus_fut);tauufut= 1*tauufut;%Qnet
% tauuhist= (hfss_hist+hfls_hist)*-1 +(rlds_hist-rlus_hist)+(rsds_hist-rsus_hist);tauuhist= -1*tauuhist;%ocean dynamics
% tauufut = (hfss_fut+hfls_fut) * -1 +(rlds_fut - rlus_fut)+(rsds_fut - rsus_fut);tauufut= -1*tauufut;%ocean dynamics

dataxx = squeeze(nanmean(tauufut,3) - nanmean(tauuhist,3));  mm=20; flag = 2;
data = nanmean(dataxx,3) ;
for i=1:size(data,1);for j=1:size(data,2);
        loc = squeeze(dataxx(i,j,:)).*data(i,j);
        loc = find(loc>=0);
        tt(i,j) = length(loc);
end;end; 

[Xq,Yq] =meshgrid(lonsst,latsst);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-200:1:200]*1,'linestyle','none'); %shading(gca,'interp')
caxis([-1 1]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:0.5:30]*mm);cbarrow

[X Y]=m_ll2xy(Xq,Yq); 
mask=tt>=30*0.90;
stipple(X,Y,mask','density',160,'color',[0 0 0]+0.5,'marker','.','markersize',5,'linewi',1);

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8); 
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

lonrange = near1(lonsst,170):near1(lonsst,360-80); latrange = near1(latsst,-5):near1(latsst,5); 
eqP   = squeeze(nanmean(nanmean(data (lonrange,latrange,:),1),2));
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);

% ------------------------------------------------------------------------
% plot regressed patterns for each terms 
close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.36 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 9.5;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

tauuhist = hfss_hist * -1;tauufut = hfss_fut *-1; %% sensitive
% tauuhist = hfls_hist * -1;tauufut = hfls_fut *-1; %% latent
% tauuhist = rlds_hist - rlus_hist;tauufut = rlds_fut - rlus_fut; % longwave
% tauuhist = rsds_hist - rsus_hist;tauufut = rsds_fut - rsus_fut; % shortwave
% tauuhist= (hfss_hist+hfls_hist)*-1 +(rlds_hist-rlus_hist)+(rsds_hist-rsus_hist);tauuhist= 1*tauuhist;%Qnet
% tauufut = (hfss_fut+hfls_fut) * -1 +(rlds_fut - rlus_fut)+(rsds_fut - rsus_fut);tauufut= 1*tauufut;%Qnet
% tauuhist= (hfss_hist+hfls_hist)*-1 +(rlds_hist-rlus_hist)+(rsds_hist-rsus_hist);tauuhist= -1*tauuhist;%ocean dynamics
% tauufut = (hfss_fut+hfls_fut) * -1 +(rlds_fut - rlus_fut)+(rsds_fut - rsus_fut);tauufut= -1*tauufut;%ocean dynamics

amoc_max_diff=[-8.51949278573794	-8.40197415486865	-5.67319697680795	-3.80951891406736	-3.96343232444736	-15.0410973871322	-13.4270047968499	-14.9252349047024	-6.56875615789509	-5.13780915944112	-10.5574300074226	-9.27216128847374	-6.18151031875871	-4.20863883172179	-7.56060870873356	-9.10073703072974	-11.0912727105611	-8.65379637434922	-6.87719377184798	-4.59023152746667	-2.31414827840000	-4.98608210086822	-4.46318017293646	-6.79151558788571	-5.02666317164445	-5.64837601777778	-15.8516325514917	-14.1840427498667	-13.6554720280000	-7.69591656083434];
im = 1:30;
for i=1:size(tauufut,1);for j=1:size(tauufut,2);
        x = (amoc_max_diff')*-1  ; 
        y = squeeze(nanmean(tauufut(i,j,:,:),3)-nanmean(tauuhist(i,j,:,:),3));
        y = y(im);
        x = x(im);
        X=[ones(length(x),1),x]; 
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;
        tt(i,j) = stats(3);     
end;end;
data = reg_moc; mm=20; flag = 3;

[Xq,Yq] =meshgrid(lonsst,latsst);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-200:1:200]*1,'linestyle','none'); %shading(gca,'interp')
caxis([-1 1]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:0.5:30]*mm);cbarrow

[X Y]=m_ll2xy(Xq,Yq); 
mask=tt<0.05;
stipple(X,Y,mask','density',160,'color',[0 0 0]+0.5,'marker','.','markersize',5,'linewi',1);

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8); 
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

lonrange = near1(lonsst,170):near1(lonsst,360-80); latrange = near1(latsst,-5):near1(latsst,5); 
eqP   = squeeze(nanmean(nanmean(data (lonrange,latrange,:),1),2));
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);


% ------------------------------------------------------------------------
% plot meridional-mean anomalies over IP basins
close all;figure(1); set(gcf, 'position', [0 100 450*1 110*1]); handaxes1 = axes('Position', [0.1 0.2 0.8 0.75]);hold on;set(gcf,'color','w');
fontsi = 9.5;
latrange = near1(latsst,-5):near1(latsst,5); 

data0 = nanmean(data(:,latrange),2)';
[h0,l0,l00]=anomaly(lonsst,smooth(data0,10),'topcolor',addcolorplus(225),'bottomcolor',addcolorplus(254),'linestyle','none');

load data_ut_thermo; load data_vt_thermo; load data_wt_thermo;
h1=plot(lonsst,smooth(data_wt_thermo,10),'LineWidth',1.2,'color',rgb('black'));
h3=plot(lonsst,smooth(data_vt_thermo,10),':','LineWidth',1.8,'color',rgb('black'));
h2=plot(lonsst,smooth(data_ut_thermo,20),'--','LineWidth',1.2,'color',rgb('black'));

ax = gca; box on;
        ax.XTick=[60:30:360];  ax.XLim =[40 290];
         ax.YTick=[-30:5:30];  ax.YLim =[-5 16];
         ax.FontSize=fontsi;
         ax.XTickLabel = {'60^oE','90^oE','120^oE','150^oE','180^o','150^oW','120^oW','90^oW'};
        ylabel(['(W m^-^2 per 8 Sv)'],'FontSize',fontsi-0.5);
        set(gca,'tickdir','out');
        set(gca,'fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.4;
        % ax.YScale='log';
ll=legend(ax,[l0 h1 h2 h3],'Q_O','Q_w_T_D','Q_u_T_D','Q_v_T_D');ll.EdgeColor ='none';ll.Color ='none';ll.FontSize=fontsi;ll.Orientation='horizontal';ll.Location='northwest';ll.Position=[0.376 0.82 0 0];
ah=axes('position',get(gca,'position'),  'visible','off');
ll2=legend(ah,[l00],'');ll2.EdgeColor ='none';ll2.Color ='none';ll2.FontSize=fontsi;ll2.Orientation='vertical';ll2.Location='southwest';ll2.Position=[0.167 0.76 0 0];











%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

varname = 'wo';
im=[];
for i=im; 
    model = modellist{i};
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    wo_hist(:,:,:,i) = nanmean(thetao,4);
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]);  thetao(thetao==0)=NaN;
    wo_fut(:,:,:,i) = nanmean(thetao,4);
end;'wo'

    latsst = ncread([filelist(1).folder,'/',filelist(1).name],'lat');
    lonsst = ncread([filelist(1).folder,'/',filelist(1).name],'lon');

rou = 1025; %kg m-3
cp = 4000; %J (K kg)-1
r0=6400000; %m
levq = [0 5,10:10:190, 200:20:280, 300:100:900, 1000:200:5800];
for j=1:size(levq,2); dz(j) = (levq(j+1)-levq(j-1))/2; end;
dx = 2*3.14*r0/360;   dy = dx;
% x - total
hist_tx = diff(to_hist,1,1)./dx; hist_tx(end+1,:,:,:) = hist_tx(end,:,:,:);
fut_tx  = diff(to_fut, 1,1)./dx; fut_tx(end+1,:,:,:)  = fut_tx(end,:,:,:);
hist_ut = hist_tx.*uo_hist;
fut_ut  = fut_tx .*uo_fut;
for i = 1:58;
hist_ut_total(:,:,i,:) = -rou*cp*hist_ut(:,:,i,:).*dz(i);
fut_ut_total(:,:,i,:) = -rou*cp*fut_ut(:,:,i,:).*dz(i);
end;
% x - due to uchange using climTx
hist_ut = hist_tx.*uo_hist;
fut_ut  = hist_tx .*uo_fut;
for i = 1:58;
hist_ut_dyn(:,:,i,:) = -rou*cp*hist_ut(:,:,i,:).*dz(i);
fut_ut_dyn(:,:,i,:) = -rou*cp*fut_ut(:,:,i,:).*dz(i);
end;
% x - due to T change using climU
hist_ut = hist_tx.*uo_hist;
fut_ut  = fut_tx .*uo_hist;
for i = 1:58;
hist_ut_thermo(:,:,i,:) = -rou*cp*hist_ut(:,:,i,:).*dz(i);
fut_ut_thermo(:,:,i,:) = -rou*cp*fut_ut(:,:,i,:).*dz(i);
end;
% y - total
hist_ty = diff(to_hist,1,2)./dy; hist_ty(:,end+1,:,:) = hist_ty(:,end,:,:);
fut_ty  = diff(to_fut, 1,2)./dy; fut_ty(:,end+1,:,:)  = fut_ty(:,end,:,:);
hist_vt = hist_ty.*vo_hist;
fut_vt  = fut_ty .*vo_fut;
for i = 1:58;
hist_vt_total(:,:,i,:) = -rou*cp*hist_vt(:,:,i,:).*dz(i);
fut_vt_total(:,:,i,:) = -rou*cp*fut_vt(:,:,i,:).*dz(i);
end;
% y - due to vchange using climTy
hist_vt = hist_ty.*vo_hist;
fut_vt  = hist_ty .*vo_fut;
for i = 1:58;
hist_vt_dyn(:,:,i,:) = -rou*cp*hist_vt(:,:,i,:).*dz(i);
fut_vt_dyn(:,:,i,:) = -rou*cp*fut_vt(:,:,i,:).*dz(i);
end;
% y - due to T change using climV
hist_vt = hist_ty.*vo_hist;
fut_vt  = fut_ty .*vo_hist;
for i = 1:58;
hist_vt_thermo(:,:,i,:) = -rou*cp*hist_vt(:,:,i,:).*dz(i);
fut_vt_thermo(:,:,i,:) = -rou*cp*fut_vt(:,:,i,:).*dz(i);
end;
% w - total
hist_tz = diff(to_hist,1,3); hist_tz(:,:,end+1,:) = hist_tz(:,:,end,:);
fut_tz  = diff(to_fut, 1,3); fut_tz(:,:,end+1,:) = fut_tz(:,:,end,:);
for i=1:58;
hist_tz(:,:,i,:) = hist_tz(:,:,i,:)./dz(i);
fut_tz(:,:,i,:) = fut_tz(:,:,i,:)./dz(i);
end;
hist_wt = hist_tz.*wo_hist;
fut_wt  = fut_tz .*wo_fut;
for i = 1:58;
hist_wt_total(:,:,i,:) = -rou*cp*hist_wt(:,:,i,:).*dz(i)*-1;
fut_wt_total(:,:,i,:) = -rou*cp*fut_wt(:,:,i,:).*dz(i)*-1;
end;
% w - due to wchange using climTz
hist_wt = hist_tz.*wo_hist;
fut_wt  = hist_tz .*wo_fut;
for i = 1:58;
hist_wt_dyn(:,:,i,:) = -rou*cp*hist_wt(:,:,i,:).*dz(i)*-1;
fut_wt_dyn(:,:,i,:) = -rou*cp*fut_wt(:,:,i,:).*dz(i)*-1;
end;
% w - due to T change using climW
hist_wt = hist_tz.*wo_hist;
fut_wt  = fut_tz .*wo_hist;
for i = 1:58;
hist_wt_thermo(:,:,i,:) = -rou*cp*hist_wt(:,:,i,:).*dz(i)*-1;
fut_wt_thermo(:,:,i,:) = -rou*cp*fut_wt(:,:,i,:).*dz(i)*-1;
end;

% ------------------------------------------------------------------------
close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.3 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 14;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

fut = fut_ut_dyn; hist = hist_ut_dyn;
% fut = fut_ut_thermo; hist = hist_ut_thermo;
% fut = fut_vt_dyn; hist = hist_vt_dyn;
% fut = fut_vt_thermo; hist = hist_vt_thermo;
% fut = fut_wt_dyn; hist = hist_wt_dyn;
% fut = fut_wt_thermo; hist = hist_wt_thermo;

im = [];
levrange = [1:5];   mm=20; % above 30m
data = nanmean(nansum(fut(:,:,levrange,im),3),4) - nanmean(nansum(hist(:,:,levrange,im),3),4); mm=20;flag =1;
dataxx = nansum(fut(:,:,levrange,im),3) - nansum(hist(:,:,levrange,im),3); 
for i=1:size(data,1);for j=1:size(data,2);
        loc = squeeze(dataxx(i,j,:)).*data(i,j);
        loc = find(loc>0);
        tt(i,j) = length(loc);
end;end; mm=20; flag = 3;

[Xq,Yq] =meshgrid(lon,lat);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-200:0.2:200]*0.2*mm,'linestyle','none'); %shading(gca,'interp')
caxis([-1 1]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:0.5:30]*mm);cbarrow

[X Y]=m_ll2xy(Xq,Yq); 
mask=tt>=30*0.9; stipple(X,Y,mask','density',160,'color',[0 0 0]+0.5,'marker','.','markersize',5,'linewi',1);

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8); 
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

lonrange = near1(lon,170):near1(lon,360-80); latrange = near1(lat,-5):near1(lat,5); 
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);

% ------------------------------------------------------------------------
close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.3 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 14;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

fut = fut_ut_dyn; hist = hist_ut_dyn;
% fut = fut_ut_thermo; hist = hist_ut_thermo;
% fut = fut_vt_dyn; hist = hist_vt_dyn;
% fut = fut_vt_thermo; hist = hist_vt_thermo;
% fut = fut_wt_dyn; hist = hist_wt_dyn;
% fut = fut_wt_thermo; hist = hist_wt_thermo;

im = [];
levrange = [1:5];   mm=20; % above 30m% 
fut  = nansum(fut(:,:,levrange,:),3);
hist = nansum(hist(:,:,levrange,:),3);
amoc_max_diff=[-8.51949278573794	-8.40197415486865	-5.67319697680795	-3.80951891406736	-3.96343232444736	-15.0410973871322	-13.4270047968499	-14.9252349047024	-6.56875615789509	-5.13780915944112	-10.5574300074226	-9.27216128847374	-6.18151031875871	-4.20863883172179	-7.56060870873356	-9.10073703072974	-11.0912727105611	-8.65379637434922	-6.87719377184798	-4.59023152746667	-2.31414827840000	-4.98608210086822	-4.46318017293646	-6.79151558788571	-5.02666317164445	-5.64837601777778	-15.8516325514917	-14.1840427498667	-13.6554720280000	-7.69591656083434];

clear tt;
for i=1:size(fut_ut,1);for j=1:size(fut_ut,2);
        x = (amoc_max_diff')*-1  ; %%%

        y = squeeze(fut(i,j,:)-hist(i,j,:));
        y = y(im);
        x = x(im);
        X=[ones(length(x),1),x]; 
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;

        tt(i,j) = stats(3);     
end;end;
data = reg_moc; flag = 3;

[Xq,Yq] =meshgrid(lon,lat);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-200:0.2:200]*0.2*mm,'linestyle','none'); %shading(gca,'interp')
caxis([-1 1]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:0.5:30]*mm);cbarrow

[X Y]=m_ll2xy(Xq,Yq); 
mask=tt<0.05;stipple(X,Y,mask','density',160,'color',[0 0 0]+0.5,'marker','.','markersize',5,'linewi',1);

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8); 
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

lonrange = near1(lon,170):near1(lon,360-80); latrange = near1(lat,-5):near1(lat,5); 
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);










