function [hout] = vocalicplot(F,varargin)

% VOCALICPLOT   Plots formant data in traditional format as a columnwise pair of plots
%
%  SYNTAX
%    H = VOCALICPLOT(F, ['Nplots', Np, 'PlotNum', ip], ['PlotHandles', hp], [...])
%    H = VOCALICPLOT(F, [...], ['MinF1', min_F1], ['MaxF1', max_F1], [...])
%
%   ARGUMENTS
%    H =  vector of axes handles
%    F =  formant matrix (expects Nx3 matrix)
%    Np = total number of plot columns in a multiplot figure (default: 1)
%    ip = index of plot column (default: 1)
%    hp = vector of plot handles for overlays (implicit hold on)
%    min_F1 = minimum plotted value for F1 (can be specified analogously for F2 and F3)
%    max_F1 = maximum plotted value for F1 (can be specified analogously for F2 and F3)
% 
%
%   DESCRIPTION
%       Plots a vertical pair of subplots, on column ip among Np columns 
%       Format is usual for vowel formants:
%       - top plot is (F1,F2), F1 vertical top-down, F2 horizontal right-to-left
%       - bottom plot is (F2,F3), F2 horizontal right-to-left, F3 vertical bottom-up
%

%
%  Created by Pierre Baraduc on 2019-07-01.
%  Copyright (C)  . All rights reserved.
%


Nplots = 1;
ip = 1;
h = [];
stdargin = {};
newplot = false;

minF1 = 200;
maxF1 = 900;
minF2 = 500;
maxF2 = 2500;
minF3 = 1500;
maxF3 = 3500;

Nvarargs = length(varargin);

k = 1;
l = 1;
while k <= Nvarargs
  stdarg = true;
  s = lower(varargin{k});
  if length(s) > 4
    switch s(1:5)
      case 'nplot'
      Nplots = varargin{k+1}; k = k+2; stdarg = false;
      case 'plotn'
      ip = varargin{k+1}; k = k+2;  stdarg = false;
      case 'ploth'
      h = varargin{k+1}; k = k+2;  stdarg = false;
      case 'minf1'
      minF1 = varargin{k+1}; k = k+2;  stdarg = false;
      case 'minf2'
      minF2 = varargin{k+1}; k = k+2;  stdarg = false;
      case 'minf3'
      minF3 = varargin{k+1}; k = k+2;  stdarg = false;
      case 'maxf1'
      maxF1 = varargin{k+1}; k = k+2;  stdarg = false;
      case 'maxf2'
      maxF2 = varargin{k+1}; k = k+2;  stdarg = false;
      case 'maxf3'
      maxF3 = varargin{k+1}; k = k+2;  stdarg = false;
    end
  end
  if stdarg
    stdargin{l} = varargin{k};
    k = k+1; 
    l = l+1;
  end
end

if ip > Nplots
  error('VOCALICPLOT: Plot number exceeds Nplots')
end

if size(F,1) == 3 
  if size(F,2) ~= 3
    F = F';
  else
    error('VOCALICPLOT: Wrong size for formant matrix (should be Nx3)')
  end
end

if isempty(h)
  newplot = true;
end
nexttile(1)
plot(F(:,2),F(:,1),stdargin{:}); hold on
if newplot 
  axis([minF2 maxF2 minF1 maxF1])
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
if ip == 1 && newplot ylabel('F1 (Hz)'), end 
end

nexttile(2)
plot(F(:,2),F(:,3),stdargin{:}); hold on
if newplot 
axis([minF2 maxF2 minF3 maxF3])
set(gca,'XDir','reverse')
xlabel('F2 (Hz)')
if ip == 1 ylabel('F3 (Hz)'), end
end

if Nplots == 1
  set(gcf,'Position',[500 200 400 700]) 
end

drawnow 

if nargout == 1 % output only when asked
  hout = h;
end

