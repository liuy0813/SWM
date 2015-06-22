%% Initialize model running parameters.

%%-------------------------------------------------------
Nens = 20;         % Number of ensemble members
Nvar = 40;         % Number of components in state variable

%%
% Observed components
% Obs = [1 6 10 16 20 26 32 40];
Obs = 1 : Nvar; % Observe all


%%
% Observe odd number of states among the first 20 states and all states from 20-40 states.
% Obs  = [1,3,5,7,9,11,13,15,17,19,21,22,23,24,25,26,27,...
%    28,29,30,31,32,33,34,35,36,37,38,39,40];

% Number of observed component
m = length(Obs);

% Observation mapping operator, H is a fat short matrix
H = zeros(m, Nvar);
for i = 1 : m
    H(i, Obs(i)) = 1;
end

%% Initial spin-up run

spinup = 10;
tspan1_spinup = 0 : 0.02 : spinup;
Anodes = 1 : Nvar;
ini = linspace(-2, 2, Nvar);  % initial condition
UnObs = setdiff(Anodes, Obs);   % Un-observed component

% True initial condition to start with
[tref, yref] = ode45('lorenz96', tspan1_spinup, ini);
ini_true = yref(length(tref), :)';

%% Model running time
tspan1 = 0 : 0.02 : 8;

%% Use the ini_true to generate reference solution
% Reference solution
[tReference, yReference] = ode45('lorenz96', tspan1, ini_true);

%% Construct background covariance matrix at t0

L = 1;  % Correlation distance
Bck_cov_t0 = zeros(Nvar, Nvar);

% use 2% perturbation as standard deviation to generate background
% covariance matrix.

pert(1 : Nvar) = 0.02 * ini_true(1 : Nvar);

for i = 1 : Nvar
    for j = 1 : Nvar
        d = abs(i - j);
        if (d > Nvar / 2)
            d = Nvar - d;
        end
        Bck_cov_t0(i, j) = pert(i) .* pert(j) .* exp( - (d / L)^2);
    end
end

%% Computing inverse of background covariance and square root of backgroud covariance

invBck_cov_t0 = inv(Bck_cov_t0);
sqBck_cov_t0 = sqrtm(Bck_cov_t0);

%% Observation points:

obs_num = length(tReference);

% ObsPoints = floor(linspace(1, length(tReference), 30));

ObsPoints = floor(linspace(1, length(tReference), obs_num));
ObsTimes  = tReference(ObsPoints); % time
ObsValues = yReference(ObsPoints, :) * H'; % mapped to the state space.

%% Start with a perturbed initial guess, 10% perturbation
ini_pert = ini_true * 1.1;

% Generate initial state ensemble using Gaussian random perturbations
ran_err = randn(Nvar, Nens);


%% A is the matrix holding state ensembles
for j = 1 : Nens
    A(:, j) = ini_pert + ran_err(:, j);   % A: Nvar-by-Nens
end


% Mean of the initial ensemble.
IniMean = mean(A, 2);

% The predicted solution
[tPredict, yPredict] = ode45('lorenz96', tspan1, ini_pert);

% The Analysis solution
tAnalysis = tspan1(1);        
yAnalysis = ini_pert';
yObservation = yAnalysis;      % Observation is generated from reference solution

