clear all;
global goalseq goalcodeseq minInput maxInput svdmatrix mu
global audioGoalMean proprioGoalMean audioGoalInvCov proprioGoalInvCov

global WE WD  bE bD xmin xmax
global Nx Nu Ns

%%%%%%%%%%%%%%%%%%%%%
java = false;
%%%%%%%%%%%%%%%%%%%%%

if java
  Ncolors = 8;
  co=colormamp;
  plotcolor = [1 0 0;
    0 1 0
    0 0 1
    0 1 1
    1 0 1
    1 1 0
    0 0 0
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560

    0.4660 0.6740 0.1880];
end

load('enc5')
WE = autoenc.EncoderWeights;
WD = autoenc.DecoderWeights;
bE = autoenc.EncoderBiases;
bD = autoenc.DecoderBiases;
xtemp = autoenc.network.inputs{1}.range;
xmin = xtemp(:,1);
xmax = xtemp(:,2);
%%

load('vocaltract_data')
load('svd_mu_simuls4.mat');

JAW_CLOSURE = 2.5; % mm jaw closure in C++ code
SELPT = 12; % mid tongue

itgt = 1;
D = 2; % delays for vel
Nx = 5; % == DOF
Nxx = 16;
Nu = 7;
Nf = 3;
Np = 3; % proprioception
Ns = 1+(D+2)*Nx;
N_SLICES = 31;
Nt = N_SLICES; % tactile FB
Ny = 43; % 3*Np + Nf + Nt; % complete feedback

load('Sa_Regions_test.mat')
load('Sa_Regions_proprio2.mat')

audioGoalMean = Mean;
proprioGoalMean = M_proprio;
audioGoalInvCov = iCov;
proprioGoalInvCov = iCov_proprio;
goalOrder = 'ie3a@Okt'; %'iemaxckt';

DT = 0.002;%

%%%%%%%%

load('System_autoenc_LSTM_bis8')
minInput = PX.xmin;
maxInput = PX.xmax;


% BASELINE MOVE %%%%%%%%%%%
%optimize_move = true;
optimize_move = false; % to save computation time if only additional noisy simulations are needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
Noptim_iter = 20000;
Nquickoptim_iter = 2500;
%%%%%%%%%%%%%

p0 = encode(autoenc, zeros(32,1));
u0 = [50,44,51,95,18,90,63];
s0 = [0; p0; zeros(Ns-Nx-1,1)];

sim_f = [];
T = 0.12;

% noise levels across conditions
unoise = [1e5*ones(1,3) 1e-6];
anoise = [1e8 1e-4 1e8  1e-6]; 
snoise = [1e8 1e8  1e-4 1e-6];
Nsimul = [20*ones(1,3) 1]; % one simulation is enough without noise


%%%%%% CHOICE OF NOISE CONDITION
for c = 1:4

  noise.var_u_add = unoise(c);
  noise.var_s_add = [anoise(c) snoise(c) snoise(c)]; %[1e-4 1e-6 1e-6];%
  noise.var_x_add = 0;
  noise.var_u_sdn = 0;
  noise.var_s_sdn = [0 0 0];
  DELAY = [0.08 0.03 0.03];

  goalcodeseq = 'ttee';
  goaltimes = [T T+0.07 T+0.12 T+0.16]; 
  Ng = length(goaltimes);
  goalmodalities = [1 1 1 1;1 1 1 1;1 1 0 0];
  dummy = zeros(Nf+Np,Ng);
  goalseq = [dummy;goaltimes; goalmodalities];


  zeroperturb.audiovector = [0 0 0];
  zeroperturb.t_onset = 0.0;
  zeroperturb.t_full = 0;



  sequence_PP = PPcoding(goalcodeseq)
  L_sequence = length(sequence_PP);
  u_goalinit = zeros(L_sequence+1,Nu);
  u_goalinit(1,:) = u0;
  for g = 1:L_sequence
    dicoV=['M_Dico[', sequence_PP(g), '].mat'];
    load(dicoV);
    u_goalinit(g+1,:) = M(1,:);
  end

  t = 0:DT:goaltimes(end);
  N = length(t);
  dt = goaltimes(end)/(N-1);

  uinit = interp1([0 goaltimes(1:end-1)-0.02 goaltimes(end)]',[u_goalinit], t)';

  init_KalmanHS = zeros(1,20);
  init_KalmanCS = zeros(1,20);
  init_PlantHS = zeros(1,20);
  init_PlantCS = zeros(1,20);
  E0 = zeros(Ns,Ns);

  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % vanilla whole-trajectory optimization
  if optimize_move
    tic
    [u_ref, s, fval] = direct_optim(@lstmmultisens_via_dop0funs_nosmoothing,uinit,s0,t,0,goaltimes,init_KalmanHS,init_KalmanCS,Noptim_iter);
    %  [u_ref, s, fval] = direct_optim(@lstmmultisens_via_dop0funs_nosmoothing,uinit,s0,t,0,goaltimes,init_KalmanHS,init_KalmanCS,Noptim_iter);
    toc
    [sim_s_ref, sim_sest_ref, sim_y_ref, sim_ynoisy_ref, sim_yest_ref, sim_noisy_u_ref, E_ref] = lstmmultisens_simnoisy_nosmoothing(u_ref', s0, s0, E0, t, zeroperturb, DELAY, noise, true, true, init_KalmanHS, init_PlantCS, init_PlantHS, init_PlantCS);
    sim_s = sim_s_ref(1:end-1,:);
    sim_sest = sim_sest_ref(1:end-1,:);
    sim_y = sim_y_ref;
    sim_ynoisy = sim_ynoisy_ref;
    save('simul_te_base', 'sim_s','sim_sest','sim_ynoisy','sim_y', 'u_ref', 'E_ref');
  else
    load('simul_te_base', 'u_ref');
    [sim_s_ref, sim_sest_ref, sim_y_ref, sim_ynoisy_ref, sim_yest_ref, sim_noisy_u_ref, E_ref] = lstmmultisens_simnoisy_nosmoothing(u_ref', s0, s0, E0, t, zeroperturb, DELAY, noise, true, true, init_KalmanHS, init_PlantCS, init_PlantHS, init_PlantCS);
  end
  %%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pseudo-continuous feedback control optimization
  Nsim = Nsimul(c);

  TIME_SLICE = 10; % time steps (= 20 ms)
  new_t = t(TIME_SLICE+1:end);
  new_N = size(new_t,2);
  Nslices = ceil(new_N/TIME_SLICE);
  res_goalseq = goalseq;
  ref_goalcodeseq = goalcodeseq;


  for nbr = 1:Nsim
    uref = u_ref;
    E0 = E_ref;
    goalseq = res_goalseq;
    goalcodeseq = ref_goalcodeseq;
    sim_s = nan(N+1,Ns);
    sim_sest = nan(N+1,Ns);
    sim_y = nan(N,Ny);
    sim_ynoisy = nan(N,Ny);
    sim_yest = nan(N,Ny);
    sim_u = nan(N,Nu);

    HS_Kalman = zeros(1,20);
    CS_Kalman = zeros(1,20);
    HS_Plant = zeros(1,20);
    CS_Plant = zeros(1,20);
    s_init = s0;
    sest_init = s0;
    slice_indx = [1:TIME_SLICE ];
    u_slice = uref(:,1:TIME_SLICE )';
    [s_slice, sest_slice, y_slice, ynoisy_slice,  yest, u_slice, Eest, HS_Kalman, CS_Kalman] = lstmmultisens_simnoisy_nosmoothing(u_slice, s_init, sest_init, E0, t(slice_indx), zeroperturb, DELAY, noise, true, true, HS_Kalman, CS_Kalman, HS_Plant, CS_Plant);
    sim_s(1:TIME_SLICE ,:) = s_slice(1:TIME_SLICE,:);
    sim_sest(1:TIME_SLICE ,:) = sest_slice(1:TIME_SLICE,:);
    sim_y(1:TIME_SLICE ,:) = y_slice(1:TIME_SLICE,:);
    sim_ynoisy(1:TIME_SLICE ,:) = ynoisy_slice(1:TIME_SLICE,:);
    sim_yest(1:TIME_SLICE ,:) = yest(1:TIME_SLICE,:);
    sim_u(1:TIME_SLICE ,:) = u_slice(1:TIME_SLICE,:);
    s_init = s_slice(end,:);
    sest_init = sest_slice(end,:);
    for k = 2:Nslices-1
      slice_indx = (k-1)*TIME_SLICE + [1:(TIME_SLICE)];
      while (t(slice_indx(1)) >= goalseq(Nf+Np+1,1))
        goalseq(:,1) = [];
        goalcodeseq(:,1) = [];
      end
      u_opt = direct_optim(@lstmmultisens_via_dop0funs_nosmoothing,uref(:,slice_indx(1):end),sest_init,t(slice_indx(1):end),1,goalseq(Nf+Np+1,:), HS_Kalman, CS_Kalman, Nquickoptim_iter); % first optimization with Nnodes = number of goals
      u_slice = u_opt(:,1:TIME_SLICE)';
      [s_slice, sest_slice, y_slice, ynoisy_slice, yest, u_slice, Eest, HS_Kalman, CS_Kalman] = lstmmultisens_simnoisy_nosmoothing(u_slice, s_init, sest_init, Eest, t(slice_indx), zeroperturb, DELAY, noise, false, false, HS_Kalman, CS_Kalman, HS_Plant, CS_Plant);
      s_init = s_slice(end,:);
      sest_init = sest_slice(end,:);
      sim_s(slice_indx,:) = s_slice(1:TIME_SLICE,:);
      sim_sest(slice_indx,:) = sest_slice(1:TIME_SLICE,:);
      sim_y(slice_indx,:) = y_slice;
      sim_ynoisy(slice_indx,:) = ynoisy_slice;
      sim_yest(slice_indx,:) = yest;
      sim_u(slice_indx,:) = u_slice;
      uref(:,slice_indx(1):end) = u_opt;
    end
    slice_indx = ((Nslices-1)*TIME_SLICE+1):length(t);
    if (slice_indx == N)
    else
      n = length(slice_indx);
      while (t(slice_indx(1)) >= goalseq(Nf+Np+1,1))
        goalseq(:,1) = [];
        goalcodeseq(:,1) = [];
      end
      [u_opt, s, fval] = direct_optim(@lstmmultisens_via_dop0funs_nosmoothing,uref(:,slice_indx(1):end),sest_init,t(slice_indx(1):end),1,goalseq(Nf+Np+1,:), HS_Kalman, CS_Kalman,Nquickoptim_iter); % first optimization with Nnodes = number of goals
      u_slice = u_opt(:,1:n)';
      [s_slice, sest_slice, y_slice, ynoisy_slice, yest, u_slice, Eest, HS_Kalman, CS_Kalman] = lstmmultisens_simnoisy_nosmoothing(u_slice, s_init, sest_init, Eest, t(slice_indx), zeroperturb, DELAY, noise, false, false, HS_Kalman, CS_Kalman, HS_Plant, CS_Plant);
      sim_s(slice_indx,:) = s_slice(1:n,:);
      sim_sest(slice_indx,:) = sest_slice(1:n,:);
      sim_y(slice_indx,:) = y_slice;
      sim_ynoisy(slice_indx,:) = ynoisy_slice;
      sim_yest(slice_indx,:) = yest;
      sim_u(slice_indx,:) = u_slice;
      uref(:,slice_indx(1):end) = u_opt;
    end

    N = size(sim_s,1)-1;
    sim_x = decode(autoenc,sim_s(:,2:Nx+1)')' + repmat(mu',N+1,1);
    sim_xest = decode(autoenc,sim_sest(:,2:Nx+1)')' + repmat(mu',N+1,1);

    %saving trajectory variables
    save(['simul_te_cond',num2str(c),'_trial',num2str(nbr)],'sim_s','sim_sest','sim_y','sim_ynoisy','sim_yest','sim_u','Nsim','u_ref');
  end

end %condition

function PPcode = PPcoding(sampacode)
PPcode = sampacode;
for i = 1:length(sampacode)
  switch sampacode(i)
    case '@', PPcode(i) = 'x';
    case 'O', PPcode(i) = 'c';
    case '3', PPcode(i) = 'm';
  end
end
end