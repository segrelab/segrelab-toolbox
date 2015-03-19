function solution = fba(model,parameters)
%FBA Flux balance analysis.
%
% Runs flux balance analysis (FBA) on the provided model.
% 
% Usage
% -----
%
% SOLUTION = FBA ( MODEL , PARAMETERS )
%
% Inputs
% ------
%
% MODEL = A flux balance model.
%
% Outputs
% -------
%
% SOLUTION = Structure. A solution structure prociding eh following fields:
%
%   ObjectiveValue
%
%   Fluxes
%
% References
% ----------
%
% 

	% TODO: validate inputs
	
	% convert model to Gurobi LP problem
	problem.A = model.S;
	problem.obj = [model.Reactions.ObjectiveCoefficient];
	problem.lb = [model.Reactions.LowerBound];
	problem.ub = [model.Reactions.UpperBound];
	problem.rhs = zeros(1,size(problem.A,1));
	problem.vtype = 'C';
	problem.sense = '=';
	problem.modelsense = 'max';
	
	% TODO: parameters
	
	% solve and create solution structure
	solution.Solver = ...
		gurobi(problem,segrelab.fba.solver.gurobi.DEFAULT_FBA_PARAMETERS);
	solution.ObjectiveValue = solution.Solver.objval;
	solution.Fluxes = solution.Solver.x;
	solution.ShadowPrices = solution.Solver.pi;
	solution.ReducedCosts = solution.Solver.rc;
	solution.IsOptimal = strcmp(solution.Solver.status,'OPTIMAL');
	
end
