load('simul_te')
load('enc5')
%% TANGENTIAL VELOCITY
t1 = t;%[0:DT:T];
data1=[];
data2=[];
N = 135;
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    sim_x = decode(autoenc,sim_s(:,2:Nx+1)')' + repmat(mu',N+1,1);
    y = squeeze(sim_x(1:end,:));
    sim_v4 = diff(y',[],2)./DT;
    vtemp1 = sqrt(sim_v4(12,:).^2 + sim_v4(28,:).^2);
    vtemp2 = sqrt(sim_v4(6,:).^2 + sim_v4(22,:).^2);
    data1 = cat(1,data1,vtemp1);
    data2 = cat(1,data2,vtemp2);
end
va1 = median(data1);
sda11 = prctile(data1, 25);
sda12 = prctile(data1, 75);
vb1 = median(data2);
sdb11 = prctile(data2, 25); 
sdb12 = prctile(data2, 75);

figure();
tiledlayout(2,1,"TileSpacing","tight")
nexttile
plot(t(1:end-1)*1000, va1,'g','LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[sda12, fliplr(sda11)],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
xlabel({'Time (ms)'},'FontSize', 22)
ylabel('Tangential velocity (mm/s)')
set(gca,'FontSize', 22);
text(10,135,'Apical region','FontSize', 24);

ylim([0 200])
nexttile
plot(t(1:end-1)*1000, vb1,'g','LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[sdb12, fliplr(sdb11)],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
xlabel({'Time (ms)'},'FontSize', 22)
ylabel('Tangential velocity (mm/s)')
set(gca,'FontSize', 22);
text(10,170,'Pharyngeal region','FontSize', 24);


%% DISPLACEMENT IN X/Y (AXIS)
Nx = 5;
load('svd_mu_simuls4')
% t(end) = [];
data1 = [];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    essai = sim_s(1:end-1,2:Nx+1)';
    dec_essai = decode(autoenc,essai);
    data1 = cat(3,data1,dec_essai);
end
data11 = median(data1(6,:,:),3);
data12 = median(data1(12,:,:),3);
data13 = median(data1(22,:,:),3);
data14 = median(data1(28,:,:),3);
sda11 = sqrt(var(squeeze(data1(6,:,:))'));
sda12 = sqrt(var(squeeze(data1(12,:,:))'));
sda13 = sqrt(var(squeeze(data1(22,:,:))'));
sda14 = sqrt(var(squeeze(data1(28,:,:))'));

figure();
tiledlayout(2,1,"TileSpacing","compact")
nexttile
plot(t(1:end-1)*1000, data11,'g','LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[prctile(squeeze(data1(6,:,:))',75), fliplr(prctile(squeeze(data1(6,:,:))',25))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
title('Pharyngeal region')
set(gca,'FontSize',24)
nexttile
plot(t(1:end-1)*1000, data13,'g','LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[prctile(squeeze(data1(22,:,:))',75), fliplr(prctile(squeeze(data1(22,:,:))',25))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
xlabel({'Time (ms)'})
set(gca,'FontSize',24)

figure()
tiledlayout(2,1,"TileSpacing","compact")
nexttile
plot(t(1:end-1)*1000, data12,'g','LineWidth',1); hold on
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[prctile(squeeze(data1(12,:,:))',75), fliplr(prctile(squeeze(data1(12,:,:))',25))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylabel('displacement in X(mm)')
title('Apical region')
set(gca,'FontSize',24)
nexttile
plot(t(1:end-1)*1000, data14,'g','LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[prctile(squeeze(data1(28,:,:))',75), fliplr(prctile(squeeze(data1(28,:,:))',25))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
xlabel({'Time (ms)'})
ylabel('displacement in Y(mm)')
set(gca,'FontSize',24)


%% FORMANTS
figure(); tiledlayout(3,1,"TileSpacing","compact")
nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = bark2hz(sim_y(1:end,1))';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));

plot(t*1000, (median1b),'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[median1b-sd1, fliplr(median1b+sd1)],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
xlabel({'Time (ms)'})
ylabel(' F1(Hz)')
set(gca,'FontSize',20)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = bark2hz(sim_y(1:end,2))';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));

plot(t*1000, median1b,'g','LineWidth',1); hold on
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[median1b-sd1, fliplr(median1b+sd1)],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;

xlabel({'Time (ms)'})
ylabel(' F2(Hz)')
set(gca,'FontSize',20)


nexttile
data1=[];
data2=[];
data3=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = (bark2hz(sim_y(1:end,3))');
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));

plot(t*1000,  (median1b),'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[median1-sd1, fliplr(median1+sd1)],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;

xlabel({'Time (ms)'})
ylabel(' F3(Hz)')
set(gca,'FontSize',20)

%% LAMBDA VALUES
figure();
tiledlayout(2,4,"TileSpacing","compact")
nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,1)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])
xlabel({'Time (ms)'})
ylabel('Lambda values (mm)')
set(gca,'FontSize',20)
legend({'GGP'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,2)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])

xlabel({'Time (ms)'})
set(gca,'FontSize',20)
legend({'GGA'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,3)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])

xlabel({'Time (ms)'})
set(gca,'FontSize',20)
legend({'HYO'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,4)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])
xlabel({'Time (ms)'})
set(gca,'FontSize',20)
legend({'STYLO'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,5)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])
xlabel({'Time (ms)'})
ylabel('Lambda values (mm)')
set(gca,'FontSize',20)
legend({'VERT'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,6)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])
xlabel({'Time (ms)'})
set(gca,'FontSize',20)
legend({'SL'},'FontSize',18)


nexttile
data1=[];
for i = 1:20
    load(['simul_te',num2str(i),'.mat'])
    y2 = sim_u(1:end,7)';
    data1 = cat(1,data1,y2(:,:));
end
median1 = median(data1);
sd1 = sqrt(var(data1));
plot(t*1000, median1,'g','LineWidth',1); hold on 
fill([t(1:end)*1000, fliplr(t(1:end)*1000)],[(prctile(data1,75)), fliplr((prctile(data1,25)))],'k','FaceColor',[.1 .9 .5],'EdgeColor','none','FaceAlpha',.2); hold on;
ylim([-20 120])
xlabel({'Time (ms)'})
set(gca,'FontSize',20)
legend({'IL'},'FontSize',18)

