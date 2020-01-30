% 14/07/2019
%inputs: number of genes, genome ID, parameters (sig, asym, adj), ID of initial condition, filename where the graph is stored
%outputs: r1(binary vector of 'pluripotent' cells that belong to the graph), frac_r1 (fraction of 'pluripotent' cells) numsteps1 (number of transition steps taken by each graph-cell-type to reach its final graph), r2(binary vector of non-graph cells that map to the input graph), frac_r2(fraction of non-graph cells that map to the input graph),numsteps2(number of steps taken by each non-graph cell to reach its final graph)
   
function [r1,frac_r1,numsteps1,r2,frac_r2,numsteps2] = explore_reproducibility(N,genomeID,initID,Sig,asym1,adj1,fname)


load('genome_initialcondition_seeds.mat');
allgenomes = A% allgenomes is a column vector containing all the seeds

load('Adj_asym_sig_seeds.mat');% A2 is a column vector containing seeds for generating rules matrices

fname1 = sprintf('%s/lineage_maps_genome%d_N%d.mat', fname, genomeID, N);
load(fname1);

if initID <= size(lineage_maps,1)
	linmap1 = lineage_maps{initID,Sig,asym1,adj1};
else
	linmap1 = [];
end

if any(linmap1)

	% initial conditions tested here are all single steady states. not testing combinations.
	numinit = 10;% number of initial conditions to be tried for each genome

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

	  r1 = false(size(s_final));
	  numsteps1 = zeros(size(s_final));

		for i1 = 1 : length(s_final)
		  disp('initialized..');
		  s1_id = s_final(i1);
		  organism_hash = zeros(1,numss);% keeps track of the cells composing the organism at all time steps
		  organism_hash(s1_id) = 1; % unique hash for each collection of cell-types
		  state_repeat = 0; % if the set of cell-types has been seen before, state_repeat=1, else =0

		  count1 = 0; % stop if code does not converge within 1000 steps  

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
		    %disp(s1_id);%CHECKING OVERLAP OF TRAJECTORY
		    s1_id = s4_id;
		    count1 = count1 + 1;

		  end
		  numsteps1(i1) = count1;
		  if isequal(s1_id,s_final)
		    r1(i1) = 1;
		  end
		end
	  frac_r1 = sum(r1)/length(r1);
	  
	% r2 = check_reproducibility(numss,B,sig,Daughters,mother_index,adjacencies,s_final,N);
	  if any(s_left)
	    disp(s_final);
	    disp(s_left);
	    disp(numss);
		  r2 = false(size(s_left));
		  numsteps2 = zeros(size(s_left));
			for i1 = 1 : length(s_left)
			  disp('initialized..');
			  s1_id = s_left(i1);
			  organism_hash = zeros(1,numss);% keeps track of the cells composing the organism at all time steps
			  organism_hash(s1_id) = 1; % unique hash for each collection of cell-types
			  state_repeat = 0; % if the set of cell-types has been seen before, state_repeat=1, else =0

			  count1 = 0; % stop if code does not converge within 1000 steps  

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
			    %disp(s1_id);%CHECKING OVERLAP OF TRAJECTORY
			    s1_id = s4_id;
			    count1 = count1 + 1;

			  end
			  numsteps2(i1) = count1;
			  if isequal(s1_id,s_final)
			    r2(i1) = 1;
			  end
			end
	  frac_r2 = sum(r2)/length(r2);
	  else 
	    r2 = 2;
	    numsteps2 = 1001;
	    frac_r2 = 0;
	  end


else
	r1 = 2;
	r2 = 2;
	numsteps1 = 1001;
	numsteps2 = 1001;
	frac_r1 = 0;
	frac_r2 = 0;
end
