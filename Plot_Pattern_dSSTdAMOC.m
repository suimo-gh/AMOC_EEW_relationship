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
        temp = ncread([filelist(1).folder,'/',filelist(1).name],'ts');
    ssthist(:,:,:,i) = temp(:,:,:);
end;clear temp;ssthist=squeeze(nanmean(ssthist,3));
for i=1:30; i
    model = modellist{i};
        filelist = dir([pathi2,'/',varname,'_',model,'_',expid2,'_run1_1x1_ann.nc']); %in 'kg s-1'
        temp = ncread([filelist(1).folder,'/',filelist(1).name],'ts');
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

varname = 'wo';
im = [ ];% wo models
for i=im; i 
    model = modellist{i};    
        filelist = dir([pathi1,'/',varname,'_',model,'_',expid1,'_run1_1x1_ann.nc']); %in 'kg s-1'
        load([filelist(1).folder,'/',filelist(1).name]); 
    wo_hist(:,:,:,i) = nanmean(thetao,4);
end;'wo'


% ------------------------------------------------------------------------
% plot dSST pattern
close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.36 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 9.5;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:4:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,ones(0,3),colorfinal(18:end,:));colormap(colorfinal);

im=1:30;
dataxx = squeeze((sstfut(:,:,im)) - (ssthist(:,:,im)));  mm=1;
data = nanmean(dataxx,3) ;

[Xq,Yq] =meshgrid(lonsst,latsst);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-200:1:200]*0.1*mm,'linestyle','none'); %shading(gca,'interp')
caxis([1.5 3]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-0.5:0.5:6]);cbarrow

mm=10^-6;
levrange = near1(lev,30);
dataw = nanmean(nanmean(wo_hist(:,:,levrange,:),3),4); % m/s
hold on;    [C1,h1]=m_contour(Xq,Yq,dataw',[0.5:1.5:1000]*3*mm,'linestyle','-','color',rgb('bright blue'),'linewi',0.5); %shading(gca,'interp');[-100:0.5:-0.5]*10^-6,

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8); 
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');set(gca,'xcolor',[0 0 0])

lonrange = near1(lonsst,170):near1(lonsst,360-90); latrange = near1(latsst,-5):near1(latsst,5); 
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);


% ------------------------------------------------------------------------
% plot dAMOC pattern
im = 1:30;
MOChist= moc_at_hist(:,:,im); MOCfut = moc_at_fut(:,:,im); flag = 1;
data = nanmean(MOCfut,3) - nanmean(MOChist,3); mm=1; 

close all;figure(1); set(gcf, 'position', [0 100 450*1 160*1]); handaxes1 = axes('Position', [0.1 0.36 0.8 0.6]);hold on;set(gcf,'color','w');
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);
fontsi = 9.5;

[latp,pfullp] =meshgrid(latq,levq);
hold on;    [C0,h0]=contourf(latp,pfullp,data',[-20:0.5:20]*mm,'linestyle','none','linecolor',rgb('white'),'linewidth',0.1);  
caxis([-2 2]*5*mm); 
c=colorbar('horizontal','position',[0.1 0.16 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-30:5:30]*mm);cbarrow

clim = nanmean(MOChist,3);
hold on;    [C1,h1]=contour(latp,pfullp,clim',[3:3:30],'linestyle','-','linecolor',rgb('dark grey'),'linewidth',0.8);  
clabel(C1,h1,[6:6:18],'color',rgb('black'),'LabelSpacing',400)

ax = gca; box on;
        ax.XTick=[-30:10:70]; ax.XLim =[-30 70];ax.XTickLabel={'30^oS','20^oS','10^oS','0^o','10^oN','20^oN','30^oN','40^oN','50^oN','60^oN','70^oN'};
        ax.FontSize=fontsi;
         ax.YTickLabel = {''};
         ylabel('Depth (m)','FontSize',fontsi-0.5);
        ax.YTickLabel={'10','1000' '2000' '3000' '4000','5000'}; 
        ax.YTick=[10 1000 2000 3000 4000 5000];
        ax.YLim =[10 4000]
        set(gca,'tickdir','out');
        set(gca,'ydir','reverse','fontsize',fontsi,'linewidth',1,'xminortick','off','yminortick','off');%'YScale','log'       
        ax.Box='on';ax.LineWidth=1;ax.TickLength=[0.02 0.02]*0.4;

% ------------------------------------------------------------------------
% plot regressd dSST 
im = 1:30; 
    latrange = near1(latq,-00):near1(latq,65); 
    levrange = near1(levq,500):near1(levq,2000); clear amoc;
    for i=1:size(moc_at_fut,3); 
        amoc_max_hist(i) = max(max(moc_at_hist(latrange,levrange,i)));
        amoc_max_fut(i)  = max(max(moc_at_fut (latrange,levrange,i)));
    end;
        amoc_max_diff = amoc_max_fut-amoc_max_hist;

hist= ssthist;
fut = sstfut;
for i=1:size(hist,1);for j=1:size(hist,2);
        x = (amoc_max_diff' * -1) ; %%%
        x = x(im);
        X=[ones(length(x),1),x];
        y = squeeze(fut(i,j,:)-hist(i,j,:));
        y=y(im);
        [b,bint,r,rint,stats] = regress(y,X) ;
        reg_moc (i,j)  = b(2)*8; warning off;

        tt(i,j) = stats(3);
end;end;
data =reg_moc; mm=1;

close all;figure(1); set(gcf, 'position', [0 100 450*1 120*1]); handaxes1 = axes('Position', [0.1 0.36 0.8 0.7]);hold on;set(gcf,'color','w');
m_proj('Equidistant Cylindrical','lon',[40 290],'lat',[-20 20]); fontsi = 9.5;
colorfinal = addcolorplus(275);colorfinal =colorfinal(1:5:end,:);colorfinal1 = addcolorplus(272);colorfinal1 =colorfinal1(1:5:end,:);colorfinal=cat(1,flip(colorfinal1,1),colorfinal);colorfinal = cat(1,colorfinal(2:12,:),colorfinal(16:end,:));colormap(colorfinal);

[Xq,Yq] =meshgrid(lonsst,latsst);
hold on;    [C1,h1]=m_contourf(Xq,Yq,data',[-20:0.05:20],'linestyle','none'); %shading(gca,'interp')
caxis([-0.8 0.8]*mm*1); 
c=colorbar('horizontal','position',[0.1 0.25 0.66 0.04]);  c.Label.Rotation=-90;c.Label.VerticalAlignment='bottom';c.FontSize=fontsi;c.YTick=([-4:0.4:4]*mm);cbarrow

[X Y]=m_ll2xy(Xq,Yq); 
mask=tt<=0.05;
stipple(X,Y,mask','density',160,'color',[0 0 0]+0.5,'marker','.','markersize',5,'linewi',1);

m_coast('patch',[0 0 0]+0.8,'edgecolor',[0 0 0]+0.6);  
m_coast('linewidth',1,'color',[0 0 0]+0.8);  
m_grid('box','on','linewidth',1,'linest','none','xtick',([-360:30:360]),'ytick',([-80:20:90]),'fontsize',fontsi, 'tickdir','out');

% 
lonrange = near1(lonsst,170):near1(lonsst,360-80); latrange = near1(latsst,-5):near1(latsst,5); 
kuang = zeros(size(data)); kuang (lonrange,latrange)=100;
hold on;    [C1,h1]=m_contour(Xq,Yq,kuang',[100 100],'linestyle','-','linecolor','k','linewi',0.5);



