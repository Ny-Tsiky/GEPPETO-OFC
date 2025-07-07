function format_figure_acoustics_bark(ax, formantStr,varargin)
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
    case 'F1F2'
        set (gca, 'XDir', 'reverse')
        set (gca, 'YDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [8 15], 'YLim', [1 8])
        end
        xlabel ('F2 [Bark]','FontSize',16)
        ylabel ('F1 [Bark]','FontSize',16)
    case 'F2F3'
        set (gca, 'XDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [8 15], 'YLim', [11 18])
        end
        xlabel ('F2 [Bark] ','FontSize',16)
        ylabel ('F3 [Bark]','FontSize',16)
        case 'F1F3'
        set (gca, 'YDir', 'reverse')
        if ConstrAxLim
        set (gca, 'XLim', [11 18], 'YLim', [1 8])
        end
        ylabel ('F1 [Bark]','FontSize',16)
        xlabel ('F3 [Bark]','FontSize',16)
    otherwise
        disp(['formantStr ' formantStr ' unknown ...'])
        return;
end
