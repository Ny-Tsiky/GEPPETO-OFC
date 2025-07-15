cmap = colormap('parula');
condcolor = cmap([1, 128, 256],:);
condchoice = 1:3;

Nsim = 20;

Nsurfnodes = 16;
zone_index = [6 12];

%% DATA LOADING
for c = condchoice
for i = 1:Nsim
    load(['simul_te_cond',num2str(c),'_trial',num2str(i),'.mat'])
    simx{c,i} = decode(autoenc,sim_s(1:end-1,2:Nx+1)')' + repmat(mu',N,1);
    F{c,i} = bark2hz(sim_y(1:end,1:3));
    lambda{c,i} = sim_u;
end
end

%% TANGENTIAL VELOCITY
pointName = {'Apical', 'Pharyngeal'};
figure();
tiledlayout(2,1,"TileSpacing","tight")

for c = condchoice
data{1}=[];
data{2}=[];
N = length(t);
for i = 1:Nsim
    sim_v = diff(simx{c,i}',[],2)./DT;
    for k = 1:2
      vtemp = sqrt(sim_v(zone_index(k),:).^2 + sim_v(zone_index(k)+Nsurfnodes,:).^2);
      data{k} = cat(1,data{k},vtemp);
    end
end
for k = 1:2
  med{c,k} = median(data{k});
  q1{c,k} = prctile(data{k}, 25);
  q3{c,k} = prctile(data{k}, 75);
end
end

for k = 1:2
nexttile
for c = condchoice
plot(t(1:end-1)*1000, med{c,k},'Color',condcolor(c,:),'LineWidth',1); hold on 
fill([t(1:end-1)*1000, fliplr(t(1:end-1)*1000)],[q3{c,k}, fliplr(q1{c,k})],'k','FaceColor',condcolor(c,:),'EdgeColor','none','FaceAlpha',.2); hold on;
end
axis([0 max(t)*1000 0 200])
xlabel({'Time (ms)'},'FontSize', 22)
ylabel('Tangential velocity (mm/s)')
set(gca,'FontSize', 22);
text(10,135,[pointName{k},' region'],'FontSize', 24);
end


%% DISPLACEMENT IN X/Y (AXIS)
Nx = 5;
load('enc5')
load('svd_mu_simuls4')
% t(end) = [];
for c = condchoice
  data{c} = [];
  for i = 1:Nsim
    data{c} = cat(3,data{c},simx{c,i});
  end
  for k = 1:2
    Xdata{c,k} = median(data{c}(:,zone_index(k),:),3);
    Ydata{c,k} = median(data{c}(:,zone_index(k)+Nsurfnodes,:),3);
    Xq1{c,k} = prctile(squeeze(data{c}(:,zone_index(k),:))',25);
    Yq1{c,k} = prctile(squeeze(data{c}(:,zone_index(k)+Nsurfnodes,:))',25);
    Xq3{c,k} = prctile(squeeze(data{c}(:,zone_index(k),:))',75);
    Yq3{c,k} = prctile(squeeze(data{c}(:,zone_index(k)+Nsurfnodes,:))',75);
  end
end

figure();
tiledlayout(2,2,"TileSpacing","compact")
for k = 1:2
  nexttile
  for c = condchoice
    plot(t*1000, Xdata{c,k},'g','LineWidth',1); hold on
    fill([t*1000, fliplr(t*1000)],[Xq1{c,k}, fliplr(Xq3{c,k})],'k','FaceColor',condcolor(c,:),'EdgeColor','none','FaceAlpha',.2); hold on;
    xlim([0 max(t)*1000])
    title([pointName{k},' region'])
    set(gca,'FontSize',24)
  end
end

for k = 1:2
  nexttile
  for c = condchoice
    plot(t*1000, Ydata{c,k},'g','LineWidth',1); hold on
    fill([t*1000, fliplr(t*1000)],[Yq1{c,k}, fliplr(Yq3{c,k})],'k','FaceColor',condcolor(c,:),'EdgeColor','none','FaceAlpha',.2); hold on;
    xlim([0 max(t)*1000])
    xlabel({'Time (ms)'})
    set(gca,'FontSize',24)
  end
end



%% FORMANTS
figure(); tiledlayout(3,1,"TileSpacing","compact")

for c = condchoice
  for k = 1:3
    data{c,k} = [];
    for i = 1:Nsim
      data{c,k} = cat(2,data{c,k},F{c,i}(:,k));
    end
  end
end

for k = 1:3
  nexttile
  for c = condchoice
    med{c,k} = median(data{c,k},2);
    q1{c,k} = prctile(data{c,k},25,2)';
    q3{c,k} = prctile(data{c,k},75,2)';

    plot(t*1000, med{c,k},'Color',condcolor(c,:),'LineWidth',1); hold on
    fill([t*1000, fliplr(t*1000)],[q1{c,k}, fliplr(q3{c,k})],'k','FaceColor',condcolor(c,:),'EdgeColor','none','FaceAlpha',.2); hold on;
  end
  xlim([0 max(t)*1000])
  xlabel('Time (ms)')
  ylabel(['F',num2str(k),' (Hz)'])
  set(gca,'FontSize',20)
end



%% LAMBDA VALUES
figure();
tiledlayout(2,4,"TileSpacing","compact")

for c = condchoice
  for k = 1:7
    data{c,k} = [];
    for i = 1:Nsim
      data{c,k} = cat(2,data{c,k}, lambda{c,i});
    end
  end
end

muscleName = {'GGP','GGA','HYO','STYLO','VERT','SL','IL'};

for k = 1:7
  nexttile
  for c = condchoice
    med = median(data{c,k},2);
    q1 = prctile(data{c,k},25,2)';
    q3 = prctile(data{c,k},75,2)';
    plot(t*1000, med,'Color',condcolor(c,:),'LineWidth',1); hold on
    fill([t*1000, fliplr(t*1000)],[q1, fliplr(q3)],'k','FaceColor',condcolor(c,:),'EdgeColor','none','FaceAlpha',.2); hold on;
  end
  xlim([0 max(t)*1000])
  ylim([-20 120])
  xlabel({'Time (ms)'})
  ylabel('Lambda values (mm)')
  set(gca,'FontSize',20)
  legend({muscleName{k}},'FontSize',18)
end

