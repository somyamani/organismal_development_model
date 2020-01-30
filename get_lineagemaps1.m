% 2/11/2018
% get lineage maps for N = {3,4,5,6,7}

% parameters: 
% 1. signaling: probability that a gene can be secreted as a signal
% 2. adjacency: probability that a given source cell can send its signals to a given target cell
% 3. asymmetry: probability that a gene that is 1 in the mother cell is 0 in its daughter cell


clc;
close all;

page_screen_output(0);
page_output_immediately(1);


load('genome_initialcondition_seeds.mat');
allgenomes = A;% allgenomes is a column matrix containing all the seeds for constructing gene regulatory network 
load('Adj_asym_sig_seeds.mat'); % A2 is a column matrix of seeds for the random number generator for each value of p_sig ,p_asym and p_adj

% initial conditions tested here are all single steady states. not testing combinations.
numinit = 10;% number of initial conditions to be tried for each genome

prob1=0:0.1:1;% all possible values of p_sig, p_asym and p_adj


N = 6;% number of genes (change for other N)

all_cellstates=dec2bin(0 : (2^N - 1),N)-'0';% all possible cell states with N node genomes

 
for genomeID = 1 : length(allgenomes) 
  rand1 = str2double(allgenomes{genomeID});
  rand('state',rand1);% set seed for random number generator
  
  % generate the genome = basins + steady states
  disp('generating genome');% debugging

  B = randi(2^N,1,2^N);% basin
  uniqB = unique(B);
  numss = length(uniqB);
  
  genome1 = {cell(length(uniqB),1), cell(length(uniqB),1)};% genome = {{steady states}, {basins}}
  for i1 = 1 : length(uniqB)
    f1 = find(B == uniqB(i1));
    B(f1) = i1;% set of cell types in the basin of steady state i1
    
    sum_ss = 0;
    while sum_ss == 0 % to ensure that there is a steady state in every basin
     ss = logical(randi(2,length(f1),1) - 1);% 0 = basin cell types, 1 = steady states
     sum_ss = sum(ss);
    end
  
    genome1{1}{i1} = all_cellstates(f1(ss),:);
    genome1{2}{i1} = all_cellstates(f1,:);
  end
  
  if numss > numinit
    initC = randperm(numss,numinit); % indices of steady states to be used as initial conditions
  else
    initC = 1 : numss;
  end
  % genome generated
  disp('genome generated');
  
  % things to store
  lineage_maps = cell(length(initC),11,11,11);% all lineage maps found for each initial condition and parameter values
  graph_properties = []; % [genome_id, initID, Sig, asym1, adj1, max_num_daughters, properties of the lineage map, reproducible or not] 
  
  for initID = 1 : length(initC)
    count = 0;
    
    for Sig=1:length(prob1);
      for asym1=1:length(prob1)
        for adj1=1:length(prob1)
          
          count = count + 1;
          
          ruleID=A2(count);% seed for mersenne twister
          rand('state',(ruleID + rand1)/2);
          
          sig=rand(1,N) <= prob1(Sig);% identities of signalling molecules% error fixed: 29/11/2018

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

          disp('daughters generated');% debugging

          cell_adj = prob1(adj1);% probability that a given cell type receives signals from another given cell type
          adjacencies = rand(2^N) <= cell_adj; % signaling adjacency matrix

          organism_hash = zeros(1,numss);% keeps track of the cells composing the organism at all time steps
          
          s1_id = initC(initID);
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

            s1_id = s4_id;
            count1 = count1 + 1;
            
          end

          per1 = isequal(organism_hash(1,:),organism_hash(2,:));

          disp('reached final organism');% debugging

          if per1 == 1
          
            % lineage map (edge list between different steady states present in the final organism)
            linmap = [];
            max_daughters = 0;% maximum number of daughters cells produced by any cell in the final organism 
            for i1 = 1 : length(s1_id)
              
              daughters = Daughters(mother_index == s1_id(i1),:); % daughters produced by the i1th celltype in the final organism
              max_daughters = max(max_daughters,size(daughters,1));
              signals_received2 = signals_received(ismember(s2,daughters,'rows'),:); % signals received by daughters of cell i1
              sig_add = (daughters + signals_received2) > 0;% state of cells right after signal exchange
              final_cells = unique(B(ismember(all_cellstates,sig_add,'rows')));% updated set of cells after one round of division, signaling and genome regulation
              edges = [s1_id(i1) * ones(length(final_cells),1), final_cells'];
              
              linmap = [linmap;edges];
            end
            
            disp('lineage map generated');% debugging

            lineage_maps(initID,Sig,asym1,adj1) = linmap;
            
            % graph properties of the lineage map
            prop1 = get_graph_prop(linmap);% 11 columns: see code get_graph_prop.m
            prop1 = [repmat([genomeID,initID,Sig,asym1,adj1,max_daughters],size(prop1,1),1),prop1];% added columns for parameter values
            % check reproducibility (i.e. pluripotency)
            rep = check_reproducibility(numss,B,sig,Daughters,mother_index,adjacencies,s1_id,N);
            prop1 = [prop1,rep * ones(size(prop1, 1),1)];
            graph_properties = [graph_properties;prop1];
            
            disp('got graph properties');% debugging

            disp ([initID,Sig,asym1,adj1]);
          end


        end
      end
    end
  end
  
% save lineage maps and graph properties
fname = sprintf('lineage_maps_genome%d_N%d.mat', genomeID, N);
save(fname,'lineage_maps','graph_properties');
disp(fname);% debugging
end


