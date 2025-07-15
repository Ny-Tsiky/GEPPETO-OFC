function [m, s, fval] = direct_optim(fun_handle,u,s0,tspan,boolperturb,goaltimes,HS0,CS0,nIter)
% direct optimization
% [m,s,fv] = direct_optim(fun_handle,u,s0,tspan,boolperturb,goaltimes,HS0,CS0,nIter)
%
% Inputs: 
% fun_handle: handle of cost function 
% u: initial command guess (Ncmds x Nsteps matrix)
% s0: initial state (Ns vector)
% tspan=time vector (length Nsteps)
% boolperturb=boolean for sensory perturbation
% goaltimes: timings for goals
% HS0, CS0: initial states of LSTM
% nIter: max #iterations
%
% Outputs
% m: final command matrix after optimization run
% s: trajectory for final command matrix
% fval: final value of cost function

% revised PB


global goal goalseq goalcodeseq
global audioGoalMean proprioGoalMean audioGoalInvCov proprioGoalInvCov
persistent fignum

if nargin < 9
    nIter = 1000;
end


ns = length(s0);
[nc, N] = size(u);
s = zeros(ns,N+1);
it = 0;
dt = 0.002; %(tspan(2)-tspan(1));
s(:,1) = s0;
% lower and upper bounds on command 
LB = [22,24,24,60,12,75,45]';
UB = [50,44,51,95,18,90,63]' + 10;
u0 = u(:,1); 


if nIter == 0 % if no max iteration number given
  options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',40000,'MaxIter',60000,'PlotFcns',@optimplotfval); 
else
 % options = optimset('TolFun',1e-4,'TolX',1e-2,'MaxFunEvals',nIter*2,'MaxIter',nIter,'PlotFcns',@optimplotfval);  
  options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',nIter*2,'MaxIter',nIter);%,'PlotFcns',@optimplotfval);  
end


fun_handle(0,0,0,0,0,0,0,0,4,0); % set goals and do some alloc


NODE_DT = [0.100 0.050 0.020 0.010];%0.05;  % choice of time intervals for parametrization; can allow progressive time refining


%  fminsearch 7*Nnodes parameters
function f = fun(x)
  persistent dummy

  m =  interp1([tspan(1) nodetimes]',[u0 reshape(x,nc,[])]', tspan)';
  if isempty(dummy), dummy = m; end
  fun_handle(0,0,0,0,HS0,CS0,dt,tspan(1),6,boolperturb); % initialize LSTM states
  out_HS = HS0; out_CS = CS0;
  s = fun_handle(dummy,m,s,s,out_HS,out_CS,dt,tspan,8,boolperturb);  % call with flag 8 takes care of all the time sequence in the C code
  f = s(1,N+1);  % integrated cost over trajectory
end

nodetimes = uniquetol(sort([tspan(1)+NODE_DT(3):NODE_DT(3):goaltimes(end), goaltimes]),1e-3); % goal times added to time parametrization if necessary
  
Nnodes = length(nodetimes);
lb = repmat(LB, 1,Nnodes);
ub = repmat(UB, 1,Nnodes);
u_optimstart = min(ub,max(u(:,round((nodetimes-tspan(1))/dt)+1),lb));

[x,fval] = fminsearchbnd(@fun,u_optimstart(:),lb(:),ub(:),options); 
m = interp1([tspan(1) nodetimes]',[u0 reshape(x,nc,[])]', tspan)';

if nargout > 1
  out_CS=CS0; % just for alloc
  out_HS=HS0;
  s = fun_handle(m,m,s,s,out_HS,out_CS,dt,tspan,8,boolperturb);
end

fun_handle(0,0,0,0,0,0,0,0,5,0); % dealloc and/or save final state

end
