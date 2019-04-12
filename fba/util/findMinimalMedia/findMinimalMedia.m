function [required_rxns] = findMinimalMedia(model,min_biomass)
%% minimal media finding algorithm
% author: Joshua Goldford
% date: 11-11-2016

% inputs: model := cobra model
%         min_biomass := a minimal growth rate
%        
% output: required_rxns := indicies of required exchange reactions

% set the min biomass to some specified value.  
% NOTE: If the objective function in the  the in the "c" field of the cobra model is not set to only the biomass reaction, 
% this will have to be commented out, and you will have to enter this into the lower bound manually.
model.lb(~~(model.c)) = min_biomass;


% find exchange reacitions, and set indicies for binary variables for MILP
U = findExcRxns(model);
idx.exchange = find(U);
idx.binary = length(model.rxns) + 1:length(model.rxns)+length(idx.exchange);
numVars.z = length(idx.binary);
numVars.v = length(model.rxns);
numVars.x = length(model.rxns)+length(idx.exchange);

% construct stoichiometric constraints
A0 = [model.S,zeros(length(model.mets),length(idx.binary))];
b0 = model.b;
csense0 = cell2mat(arrayfun(@(x) 'E',b0,'uni',0)');

% construct binary constraints
A1 = zeros(numVars.z,numVars.x);
for i = 1:numVars.z
    A1(i,idx.exchange(i)) = -1;
    A1(i,idx.binary(i)) = model.lb(idx.exchange(i));
    b1(i) = 0;
    csense1{i}= 'L';
end
csense1 = cell2mat(csense1);


% construct the MILP structure
MILP.A = [A0;A1];
MILP.b = [b0;b1'];
MILP.csense = [csense0,csense1];
MILP.lb = [model.lb;zeros(numVars.z,1)];
MILP.ub = [model.ub;ones(numVars.z,1)];
MILP.vartype = [arrayfun(@(x) 'C', 1:numVars.v,'uni',1),arrayfun(@(x) 'B', 1:numVars.z,'uni',1)];
MILP.osense = +1;
MILP.c = [zeros(numVars.v,1);ones(numVars.z,1)];

% find an initial solution via traditional FBA
x = optimizeCbModel(model);
ex = find(x.x(idx.exchange) < 0);
e = zeros(numVars.z,1);
e(ex) = 1;
x0 = [x.x;e];

% append a starting point for the MILP algorithm
MILP.x0 = x0;

% solve the cobra MILP problem
y = solveCobraMILP(MILP);

% output reaction indicies for required reactions
required_rxns = idx.exchange(~~y.int);

end


