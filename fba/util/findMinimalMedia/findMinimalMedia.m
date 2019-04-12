function [required_rxns] = findMinimalMedia(model,min_biomass,include_exchRxns,exclude_exchRxns)
%%findMinimalMedia MILP to find a minimal media set for a metabolic model
% author: Joshua Goldford
% date: 11-11-2016
%
% [required_rxns] = findMinimalMedia(model)
% [required_rxns] = findMinimalMedia(model,min_biomass,include_exchRxns,exclude_exchRxns)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.b:     Right hand side = dx/dt
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%   model.ub:    Upper bounds
%   model.rxns*:  Reaction names (*only required if using optional input
%                                 include_exchRxns*)
%
%OPTIONAL INPUTS
% min_biomass: Minimal growth rate (default=0.1)
% include_exchRxns: Exchange reactions to include by default (default=empty)
%                   name in model.rxns
% exclude_exchRxns: Exchange reactions to exclude by default (default=empty)
%                   name in model.rxns
%
%OUTPUT
% required_rxns: Indicies of required exchange reactions
%
% Meghan Thommes 04/12/2019 - Changed solver to gurobi explicitly, include
%                             or exclude specific exchange reactions

%% Check Inputs

if (nargin < 1)
    error('myfuns:findMinimalMedia:NotEnoughInputs', ...
        'Not enough inputs: need a model');
else
    if ~isstruct(model)
        error('myfuns:findMinimalMedia:IncorrectType', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S') || ~isfield(model,'b') || ~isfield(model,'c') || ~isfield(model,'lb') || ~isfield(model,'ub')
        error('myfuns:findMinimalMedia:IncorrectType', ...
            '"model" needs "S", "b", "c", "lb", and "ub" fields');
    end
end
%%

if ~exist('min_biomass','var') || isempty(min_biomass)
    min_biomass = 0.1;
end

% Gurobi Parameters
if ~exist('params','var')
    params.FeasibilityTol = 1e-9; % all values must be satisfied to this tolerance (min value)
    params.OutputFlag = 0; % silence gurobi
    params.DisplayInterval = 1; % frequency at which log lines are printed (in seconds)
else
    if ~isfield(params,'FeasibilityTol')
        params.FeasibilityTol = 1e-9; % all values must be satisfied to this tolerance (min value)
    end
    if ~isfield(params,'OutputFlag')
        params.OutputFlag = 0; % silence gurobi
    end
    if ~isfield(params,'DisplayInterval')
        params.DisplayInterval = 1; % frequency at which log lines are printed (in seconds)
    end
end

%% Set Up MILP

% set the min biomass to some specified value.  
model.lb(~~(model.c)) = min_biomass;

% find exchange reactions
idx.exchange = findExchRxns(model);

% if specifying reactions to include
if exist('include_exchRxns','var') && ~isempty(include_exchRxns)
    [~,idx.include,~] = intersect(model.rxns(idx.exchange),include_exchRxns);
    if isempty(idx.include)
        [~,idx.include,~] = intersect(model.rxnNames(idx.exchange),include_exchRxns);
        if isempty(idx.include)
            warning('myfuns:findMinimalMedia:IncorrectType', ...
                'Cannot find specified include_exchRxns');
            idx.include = [];
        end
    end
else
    idx.include = [];
end

% if specifying reactions to exclude
if exist('exclude_exchRxns','var') && ~isempty(exclude_exchRxns)
    [~,idx.exclude,~] = intersect(model.rxns(idx.exchange),exclude_exchRxns);
    if isempty(idx.exclude)
        [~,idx.exclude,~] = intersect(model.rxnNames(idx.exchange),exclude_exchRxns);
        if isempty(idx.exclude)
            warning('myfuns:findMinimalMedia:IncorrectType', ...
                'Cannot find specified exclude_exchRxns');
            idx.exclude = [];
        end
    end
else
    idx.exclude = [];
end

% update lower bounds
model.lb(idx.exchange) = -1000;

% set indicies for binary variables for MILP
idx.binary = length(model.rxns) + 1:length(model.rxns)+length(idx.exchange);

numVars.e = length(idx.binary);
numVars.m = length(model.mets);
numVars.n = length(model.rxns);
numVars.x = length(model.rxns)+length(idx.exchange);

% LP Model
LP_model.A = model.S; % stoichiometric matrix: linear constraint matrix [sparse matrix, mets x rxns]
LP_model.obj = model.c; % linear objective vector for each each col of A (rxn in S) [dense vector]
LP_model.rhs = model.b; % right-hand side vector for the linear constraints for each row of A (met in S) [dense vector]
LP_model.lb = model.lb; % lower bounds for each col of A (rxn in S) [dense vector]
LP_model.ub = model.ub; % upper bounds for each col of A (rxn in S) [dense vector]
LP_model.modelsense = 'max'; % maximize objective function
LP_model.sense = '='; % sense of the linear constraints for each row of A (met in S) [char array]
LP_model.vtype = 'C'; % continuous variables

% find an initial solution via traditional FBA
sol = gurobi(LP_model,params);
ex = find(sol.x(idx.exchange) < 0);
e = zeros(numVars.e,1);
e(ex) = 1;
x0 = [sol.x; e];

% construct binary constraints
A = zeros(numVars.e,numVars.x);
for i = 1:numVars.e
    A(i,idx.exchange(i)) = -1;
    A(i,idx.binary(i)) = model.lb(idx.exchange(i));
end

% reactions to include
lb = zeros(numVars.e,1);
if sum(idx.include)>0
    lb(idx.include) = 1;
end

% reactions to exclude
ub = ones(numVars.e,1);
if sum(idx.exclude)>0
    ub(idx.exclude) = 0;
end

% MILP Model
MILP_model.A = zeros(numVars.m+numVars.e,numVars.x); % stoichiometric matrix [mets+binary x rxns+binary]
MILP_model.A(1:numVars.m,1:numVars.n) = LP_model.A; % continuous
MILP_model.A(numVars.m+1:end,:) = A; % binary
MILP_model.A = sparse(MILP_model.A); % must be a sparse matrix
MILP_model.obj = [zeros(numVars.n,1); ones(numVars.e,1)]; % linear objective vector for each each col of A [dense vector, rxns+binary]
MILP_model.rhs = [model.b; zeros(numVars.e,1)]; % right-hand side vector for the linear constraints for each row of A [dense vector, mets+binary]
MILP_model.lb = [model.lb; lb]; % lower bounds for each col of A [dense vector, rxns+binary]
MILP_model.ub = [model.ub; ub]; % upper bounds for each col of A [dense vector, rxns+binary]
MILP_model.modelsense = 'min'; % minimize objective function
MILP_model.sense = [repmat('=',1,numVars.m), repmat('<',1,numVars.e)]; % sense of the linear constraints for each row of A [char array, rxns+binary]
MILP_model.vtype = [repmat('C',1,numVars.n), repmat('B',1,numVars.e)]; % type of variables for each row of A [char array, rxns+binary]

% append a starting point for the MILP algorithm
MILP_model.x0 = x0;

%% Solve the MILP

y = gurobi(MILP_model,params);
y.int = y.x(MILP_model.vtype == 'B');

% output reaction indicies for required reactions
required_rxns = idx.exchange(~~y.int);









end


