function [m, s, fval] = dop0_PB_direct_noisy(fun_handle,u,s0,tspan,boolperturb,goaltimes,HS0,CS0,nIter)

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
dt = 0.002; 
s(:,1) = s0;
LB = [22,24,24,60,12,75,45]';  
UB = [50,44,51,95,18,90,63]' + 10;
u0 = u(:,1);


if nIter == 0
  options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',40000,'MaxIter',60000,'PlotFcns',@optimplotfval);
else
  options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',nIter*2,'MaxIter',nIter); %,'PlotFcns',@optimplotfval);  
end

fun_handle(0,0,0,0,0,0,0,0,4,0); 
SUBGOAL_DT = [0.100 0.050 0.020 0.010];%0.05;
Npass = length(SUBGOAL_DT);

function f = fun(x)
  persistent dummy
  m =  interp1([tspan(1) subgoaltimes]',[u0 reshape(x,nc,[])]', tspan)'; 
  if isempty(dummy), dummy = m; end
  fun_handle(0,0,0,0,HS0,CS0,dt,tspan(1),6,boolperturb); 
  out_HS = HS0; out_CS = CS0;
  s = fun_handle(dummy,m,s,s,out_HS,out_CS,dt,tspan,8,boolperturb);  
  f = s(1,N+1); 
end

for pass = Npass 
  subgoaltimes = uniquetol(sort([tspan(1)+SUBGOAL_DT(3):SUBGOAL_DT(3):goaltimes(end), goaltimes]),1e-3);
  Nsubg = length(subgoaltimes);
  u_optimstart = u(:,round((subgoaltimes-tspan(1))/dt)+1)./UB;
  lb = repmat(LB, 1,Nsubg);
  ub = repmat(UB, 1,Nsubg);

  [x,fval] = fminsearchbnd(@fun,u_optimstart(:),lb(:),ub(:),options); 
  m = interp1([tspan(1) subgoaltimes]',[u0 reshape(x,nc,[])]', tspan)'; 
end

fun_handle(0,0,0,0,0,0,0,0,5,0); 
end
