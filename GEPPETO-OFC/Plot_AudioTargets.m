function Plot_AudioTargets(Plan,Vect,Cibles,R,ColOpt,LW)
PhonLab_T={'i','e','\epsilon','a','oe','c','k'};
Ell=Cibles.Mean;
Sigma=Cibles.Cov;

cmap = hsv(length(Vect));
LabPlan={'F1F2','F2F3'};
DimMod=[2,1;2 3];
nP=length(Plan);
for ip=1:nP
    subplot(1,nP,ip)
    hold on    
    format_figure_acoustics(gca, LabPlan{Plan(ip)})
    hold all
    for t=Vect
        if ColOpt
            Col=cmap(t,:);
        else
            Col='k';
        end
        [X,Y]=Ellipse_Error(Ell(t,DimMod(ip,:)),Sigma(DimMod(ip,:),DimMod(ip,:),t),R);
        h=plot(X,Y,'Color',Col,'LineWidth',LW);
        h_t=text( Ell(t,DimMod(ip,1)), Ell(t,DimMod(ip,2)) ,PhonLab_T(t),'HorizontalAlignment','center');
        set(h_t,'FontSize',15)
    end
end
end
