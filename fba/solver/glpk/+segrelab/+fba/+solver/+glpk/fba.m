function solution = fba(model,parameters)
%FBA Summary of this function goes here
%   Detailed explanation goes here

	% TODO: validate inputs
	
	% convert model to GLPK LP problem
	A = model.S;
	c = [model.Reactions.ObjectiveCoefficient];
	lb = [model.Reactions.LowerBound];
	ub = [model.Reactions.UpperBound];
	b = zeros(1,size(model.S,1));
	vartype = repmat('C',1,size(model.S,2));
	ctype = repmat('S',1,size(model.S,1));
	sense = -1;
	
	% TODO: parameters
	parameters.presol = 0;
	
	% solve & create soluton structure
	[solution.Solver.xopt,solution.Solver.fmin,solution.Solver.status,solution.Solver.extra] = ...
		glpk(c,A,b,lb,ub,ctype,vartype,sense,parameters);
	solution.ObjectiveValue = solution.Solver.fmin;
	solution.Fluxes = solution.Solver.xopt;
	solution.ShadowPrices = solution.Solver.extra.lambda;
	solution.ReducedCosts = solution.Solver.extra.redcosts;
	solution.IsOptimal = isequal(solution.Solver.status,5);

end

