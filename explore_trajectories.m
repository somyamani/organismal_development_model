% 13/08/2019
%inputs: number of genes, genome ID, parameters (sig, asym, adj), ID of initial condition, filename where the graph is stored
% process: get trajectories till steady state of ~10 random initial conditions
%outputs:1.trajectories = cell array with ([fraction of cells that match target, number of cells in current state] for each time step) for each initial condition
%2.target_reached =  binary vector to say whether the initial condition reached the target or not

function [trajectories,target_reached] = explore_trajectories(N,genomeID,initID,Sig,asym1,adj1,fname)


allgenomes = load('genome_initialcondition_seeds.txt');
allgenomes = A;% allgenomes is a cell array containing all the seeds

load('Adj_asym_sig_seeds.mat');% A2 is a column matrix with seeds for generation of rules matrices

fname1 = sprintf('%s/lineage_maps_genome%d_N%d.mat', fname, genomeID, N);
load(fname1);

if initID <= size(lineage_maps,1)
	linmap1 = lineage_maps{initID,Sig,asym1,adj1};
else
	linmap1 = [];
	disp('no linmap. try lower init');
end

numinit = min(10,2^N);% number of initial conditions to be tried for each genome
trajectories = cell(numinit,1);
target_reached = zeros(numinit,1);
if any(linmap1)

	% initial conditions tested here are all single steady states. not testing combinations.

	prob1=0:0.1:1;% all possible values of p_sig, p_asym and p_adj

	paramID = (Sig-1)*(11^2) + (asym1-1)*(11) + adj1;

	all_cellstates=dec2bin(0 : (2^N - 1),N)-'0';% all possible cell states with N node genomes



	  rand1 = str2double(allgenomes{genomeID});
	  rand('state',rand1);% set seed for random number generator

	  % generate the genome = basins + steady states
	  B = randi(2^N,1,2^N);% basin
	  uniqB = unique(B);
	  numss = length(uniqB);

	  genome1 = {cell(length(uniqB),1), cell(length(uniqB),1)};% genome = {{steady states}, {basins}}
	  S = zeros(1,2^N);% identities of cell-types in steady states

	  for i1 = 1 : length(uniqB)
	    f1 = find(B == uniqB(i1));
	    B(f1) = i1;% set of cell types in the basin of steady state i1

	    sum_ss = 0;
	    while sum_ss == 0 % to ensure that there is a steady state in every basin
	     ss = logical(randi(2,length(f1),1) - 1);% 0 = basin cell types, 1 = steady states
	     sum_ss = sum(ss);
	    end

	    S(f1(ss)) = i1;

	    genome1{1}{i1} = all_cellstates(f1(ss),:);
	    genome1{2}{i1} = all_cellstates(f1,:);
	  end
	  uniqB = unique(B);
	  numinit = min(length(uniqB),numinit);
	  trajectories = cell(numinit,1);
	  target_reached = zeros(numinit,1);
	  % genome generated

	  ruleID=A2(paramID);% seed for mersenne twister
	  rand('state',(ruleID + rand1)/2);
	  sig=rand(1,N) <= prob1(Sig);% identities of signalling molecules

	  p_asym=prob1(asym1); % probability of an 'on' gene being 'off' in one of the daughters
	  daughtermatrix = rand(2^N,N) >= p_asym; % p_asym = 0 implies daughter cells are identical to mother cells, p_asym = 1 implies daughter cells can only be [ 0 0 0 ...]. use this as a filter to get identities of daughter cells

	  mother_index = [];% identities of mother cells for each daughter cell produced 
	  Daughters = []; % set of daughters produced by all cell types present 
	  for i1 = 1 : numss
	    ss = genome1{1}{i1};
	    d1 = [];
	    for i2 = 1 : size(ss,1) % in the case of oscillations, all cell-types are considered to be present, therefore, daughters of all oscillating cells are produced at the same time. in this case, even if asym = 0, >1 daughter cells are produced
	      d1 = [d1;daughtermatrix * diag(ss(i2,:))];
	    end
	    ud = unique(d1,'rows');

	    Daughters = [Daughters;ud];
	    mother_index = [mother_index; i1 * ones(size(ud,1),1)];
	  end

	  disp('daughters generated');

	  cell_adj = prob1(adj1);% probability that a given cell type receives signals from another given cell type
	  adjacencies = rand(2^N) <= cell_adj; % signaling adjacency matrix
	  disp('adjacencies generated');


	  s_final = unique(linmap1)';% set of all cell states in the organism
	  s_left = setdiff(uniqB,s_final);% set of all cell states not in s_final
	  % pick ~10 random cell types to be initial conditions
	  init_cells = randperm(length(uniqB));
	  init_cells = uniqB(init_cells(1:numinit));
	  numsteps1 = zeros(size(init_cells));

		for i1 = 1 : numinit
		  disp('initialized..');
		  s1_id = init_cells(i1);
		  organism_hash = zeros(1,numss);% keeps track of the cells composing the organism at all time steps
		  organism_hash(s1_id) = 1; % unique hash for each collection of cell-types
		  state_repeat = 0; % if the set of cell-types has been seen before, state_repeat=1, else =0
		  Ti1 = [length(intersect(s1_id,s_final)),length(s1_id)];
		  count1 = 0; % stop if code does not converge within 1000 steps  
		  %DEBUG
		  disp('trajectory number')
		  disp(i1);
		  while (state_repeat==0) && (count1 < 1000)

		    % asymmetric division
		    s2 = unique(Daughters(ismember(mother_index,s1_id),:),'rows');% set of daughter cells produced in this time step. each row is a cell type.
		    s2_dec = ismember(all_cellstates,s2,'rows');

		    % signals produced by daughters cells
		    signals_produced = s2 * diag(sig);

		    % signaling
		    cell_adjmat_sub=adjacencies(s2_dec,s2_dec);
		    signals_received = (cell_adjmat_sub') * signals_produced; % set of signals received by each of the daughter cells produced in this time step
		    s3=(s2 + signals_received) > 0;
		    s3=unique(s3,'rows');% s3 is the state of the organism right after signal exchange

		    % all cell-types transition to a steady state. for oscillatory steady states, all cell types in the oscillation are present in the next time step

		    s4_id = unique(B(ismember(all_cellstates,s3,'rows')));% s4 is the updated organism after one round of division, signaling and genome regulation
		    s4_hash = zeros(1,numss);
		    s4_hash(s4_id) = 1;

		    organism_hash=[s4_hash;organism_hash];

		    state_repeat = any(ismember(organism_hash(1,:), organism_hash(2:end,:), 'rows'));% the state has been seen before if state_repeat = 1
		    % DEBUG
		    disp(s1_id);%CHECKING OVERLAP OF TRAJECTORY

		    s1_id = s4_id;
		    Ti1 = [Ti1;[length(intersect(s1_id,s_final)),length(s1_id)]];
		    count1 = count1 + 1;

		  end
		  target_reached(i1) = isequal(s1_id,s_final);
		  %DEBUG
		  disp('target reached : ')
		  disp(isequal(s1_id,s_final))
		  numsteps1(i1) = count1;
		  trajectories{i1} = Ti1;
		end
end
