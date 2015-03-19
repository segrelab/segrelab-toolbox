function model = sbml2model(sbml,varargin)
%SBML2MODEL Convert an SBML structure into a Segrè lab FBA model.
%
% Converts an SBML structure into an FBA model object or structure.
%
% Usage
% -----
%
% MODEL = SBML2MODEL ( SBML , [ATTRIBUTES] , ['struct'] )
%
% Inputs
% ------
%
% SBML = An SBML structure.
%
% Outputs
% -------
%
% MODEL = A flux balance model.
%
% References
% ----------
%
% See also TRANSLATESBML, SEGRE2COBRA

	% validate inputs
	%TODO
	
	% basic model information
	model.ID = sbml.id;
	model.Name = sbml.name; % TODO: name almost always empty?
	model.MIRIAM = parseAnnotation(sbml.annotation);
	
	% more complex model components: each in their own function
	compartmentNameByID = containers.Map(...
		{sbml.compartment.id,''},{sbml.compartment.name,char.empty});
	[model.Metabolites,metaboliteIndexByID] = ...
		parseSpecies(sbml.species,compartmentNameByID);
	[model.Reactions,reactionIndexByID] = ...
		parseReactions(sbml.reaction,compartmentNameByID);
	model.S = parseStoichiometries(sbml.reaction,metaboliteIndexByID);
	
	% add-on SBML package information
	if (sbml.SBML_level == 3) && isfield(sbml,'fbc_version')
		% correct reactions for flux balance constraints
		model = parseFBC(sbml,model,reactionIndexByID);
	end
	
	% turn into an object
	%TODO

end

function [metabolites,metaboliteIndexByID] = ...
		parseSpecies(species,compartmentNameByID)
%PARSESPECIES Parse SBML structure species into model metabolites.

	% skip boundary species: not included in our models
	metaboliteSpecies = species(~[species.boundaryCondition]);
	M = numel(metaboliteSpecies);
	
	% basic metabolite information
	ids = {metaboliteSpecies.id};
	names = {metaboliteSpecies.name};
	compartments = cellfun(@(x) compartmentNameByID(x),...
		{metaboliteSpecies.compartment},'UniformOutput',false);
	
	% some metabolite information "hidden" in the notes/annotation field
	[formulas] = cellfun(@(x) parseNotes(x,'FORMULA'),...
		{metaboliteSpecies.notes},'UniformOutput',false);
	[miriams] = cellfun(@parseAnnotation,...
		{metaboliteSpecies.annotation},'UniformOutput',false);
	
	% assemble information into metabolites
	metabolites = struct('ID',ids,...
		'Name',names,'Compartment',compartments,'Formula',formulas,...
		'MIRIAM',miriams);
	metaboliteIndexByID = containers.Map(ids,num2cell(1:M));

end

function [reactions,reactionIndexByID] = ...
		parseReactions(reactions,compartmentNameByID)
%PARSEREACTIONS Parse SBML structure reactions into model reactions.

	% basic reaction information
	N = numel(reactions);
	ids = {reactions.id};
	names = {reactions.name};
	reversibilities = {reactions.reversible};
	compartments = repmat({''},1,N);
	if isfield(reactions,'compartment')
		compartments = cellfun(@(x) compartmentNameByID(x),...
			{reactions.compartment},'UniformOutput',false);
	end
	
	% reaction constraints information
	[lb,ub,c,v] = cellfun(@parseKineticLaw,...
		{reactions.kineticLaw},'UniformOutput',false);
	
	% some reactions information "hidden" in the notes/annotation field
	[geneAssociations] = cellfun(@(x) parseNotes(x,'GENE_ASSOCIATION'),...
		{reactions.notes},'UniformOutput',false);
	[miriams] = cellfun(@parseAnnotation,...
		{reactions.annotation},'UniformOutput',false);
	
	% assemble information into reactions
	reactions = struct('ID',ids,...
		'Name',names,'Compartment',compartments,...
		'LowerBound',lb,'UpperBound',ub,'ObjectiveCoefficient',c,'Flux',v,...
		'Reversible',reversibilities,...
		'MIRIAM',miriams);
	reactionIndexByID = containers.Map(ids,num2cell(1:numel(reactions)));

end

function stoichiometries = parseStoichiometries(reactions,metaboliteIndexByID)
%PARSESTOICHIOMETRIES Parse SBML structure elements for stoichiometric matrix.

	N = numel(reactions);
	stoichiometries = zeros(metaboliteIndexByID.length,N); % TODO
	for n = 1:N
		reaction = reactions(n);
		for reactant = reaction.reactant
			if metaboliteIndexByID.isKey(reactant.species)
				m = metaboliteIndexByID(reactant.species);
				stoichiometries(m,n) = -reactant.stoichiometry;
			end
		end
		for product = reaction.product
			if metaboliteIndexByID.isKey(product.species)
				m = metaboliteIndexByID(product.species);
				stoichiometries(m,n) = product.stoichiometry;
			end
		end
	end
	stoichiometries = sparse(stoichiometries);

end

function model = parseFBC(sbml,model,reactionIndexByID)
%PARSEFBC Parse flux balance constraints for reaction properties.

	% objective function
	objective = strcmp(sbml.fbc_activeObjective,{sbml.fbc_objective.fbc_id});
	reactions = cellfun(@(x) reactionIndexByID(x),...
		{sbml.fbc_objective(objective).fbc_fluxObjective.fbc_reaction});
	coefficients = ...
		[sbml.fbc_objective(objective).fbc_fluxObjective.fbc_coefficient];
	model.Reactions(reactions) = arrayfun(...
		@(n,c) setfield(model.Reactions(n),'ObjectiveCoefficient',c),...
		reactions,coefficients); %#ok: setfield because it's inside arrayfun
	
	% reaction bounds
	for bound = sbml.fbc_fluxBound
		n = reactionIndexByID(bound.fbc_reaction);
		switch (bound.fbc_operation)
			case 'greaterEqual'
				model.Reactions(n).LowerBound = bound.fbc_value;
			case 'lessEqual'
				model.Reactions(n).UpperBound = bound.fbc_value;
			case 'equal'
				model.Reactions(n).LowerBound = bound.fbc_value;
				model.Reactions(n).UpperBound = bound.fbc_value;
		end
	end

end

function [lb,ub,c,v] = parseKineticLaw(law)
%PARSEKINETICLAW Parse kinetic law field for reaction properties.

	% default values
%  	lb = segrelab.fba.model.Reaction.DEFAULT_LOWER_BOUND;
	lb = -Inf;
%  	ub = segrelab.fba.model.Reaction.DEFAULT_UPPER_BOUND;
	ub = Inf;
%  	c = segrelab.fba.model.Reaction.DEFAULT_OBJECTIVE_COEFFICIENT;
	c = 0;
%  	v = segrelab.fba.model.Reaction.DEFAULT_FLUX;
	v = 0;
	if isempty(law)
		return
	end
	
	% parse the kinetic law parameters
	for parameter = law.parameter
		switch parameter.id
			case 'LOWER_BOUND'
				lb = parameter.value;
			case 'UPPER_BOUND'
				ub = parameter.value;
			case 'OBJECTIVE_COEFFICIENT'
				c = parameter.value;
			case 'FLUX_VALUE'
				v = parameter.value;
		end
	end

end

function miriam = parseAnnotation(annotation)
%PARSEANNOTATION Parse annotation field for common properties.

	miriam = containers.Map('KeyType','char','ValueType','char');
	if ~isempty(annotation)
		uris = regexp(annotation,MIRIAM_PATTERN,'names');
		miriam = containers.Map({uris.collection},{uris.entity});
	end

end

function varargout = parseNotes(notes,varargin)
%PARSENOTES Parse notes field for common properties.

	varargout = cell(size(varargin));
	for n = 1:numel(varargin)
		pattern = strcat('(?<=',varargin{n},':)(.+?)(?=<(?:html:)?/p>)');
		varargout{n} = regexp(notes,pattern,'once','match');
	end

end

function pattern = MIRIAM_PATTERN()
%MIRIAM_PATTERN Regular expression pattern for MIRIAM URIs in annotations.

	pattern = ...
		'(?<=http://identifiers.org/)(?<collection>[^/]+)/(?<entity>[^/"]+)';

end
