close all
Vect=1:7;
ColOpt='k';
LW=2;
R=2.8;
Plan=[1 2];



if Plot(1)
CiblesSa=load(sprintf('Sa_Regions'));
PosMod=[546 100 995 853;23 199 1882 612];
f1=figure('Position', PosMod(length(Plan),:));
hold on
Plot_AudioTargets(Plan,Vect,CiblesSa,R,ColOpt,LW)

PaperPos=[-1.5 1 33 12];
PaperSiz=[30 13];
set(gcf, 'Renderer', 'Painters')
                set(gcf, 'PaperPosition', PaperPos); 
                set(gcf, 'PaperSize', PaperSiz); 
saveas(f1,'Sa_Regions','pdf')
end
if Plot(2)
CiblesSo=load(sprintf('PCASo3_Regions'));
f2=figure('Position',[546 100 995 853]);
hold on
Plot_ProprioTargets('So',Vect,CiblesSo,R,ColOpt,LW)
saveas(f2,'So_Regions','png')
end
if Plot(3)
CiblesSoSh=load(sprintf('PCASoSh3_Regions'));
f3=figure('Position',[546 100 995 853]);
hold on
Plot_ProprioTargets('SoSh',Vect,CiblesSoSh,R,ColOpt,LW)
saveas(f3,'SoSh_Regions','png')
end