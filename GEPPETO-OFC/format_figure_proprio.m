function format_figure_proprio(ax, formantStr,varargin)
% formats a figure (axes limits, etc.) for formant plots (F1/F2 or F2/F3)
%
% this function should be used either to format the plot in the F1/F2-plane
% (formantStr = 'F1F2') or format the F2/F3-plane (formantStr = 'F2F3').

% written 03/2012 by RW (SPRECHart); modified 09/2012 by RW (PILIOS)

axes(ax);
grid off
hold on
%legend (tractVersion)
%title ([date])
if nargin>2
    ConstrAxLim=varargin{1};
else
    ConstrAxLim=1;
end
switch formantStr
    case 'P1P2'
        set (gca, 'XDir', 'reverse')
        set (gca, 'YDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [-40 40], 'YLim', [-50 50])
        end
        xlabel ('PC2','FontSize',20)
        ylabel ('PC1','FontSize',20)
    case 'P2P3'
        set (gca, 'XDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [-40 40], 'YLim', [-10 15])
        end
        xlabel ('PC2','FontSize',20)
        ylabel ('PC3','FontSize',20)
        case 'P1P3'
        set (gca, 'YDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [-10 15], 'YLim', [-50 50])
        end
        ylabel ('PC1','FontSize',20)
        xlabel ('PC3','FontSize',20)
    otherwise
        disp(['formantStr ' formantStr ' unknown ...'])
        return;
end
