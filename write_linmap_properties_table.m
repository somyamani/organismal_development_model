% 16/1/2019

%output: [Number of genes,Number of fixedpoints,Number of oscillators,basin_vec,ss_vec,sig,sig_vec,adj,adjmat,asym,divmat,linmap_edgelist]

% write data into a text file. delimiters= ';' within a
% single object AND ',' between objects


N=3; number of genes

filename = 'allgraphs_N3_1.txt';
fid = fopen(filename,'a');

disp('opened file for writing');


load('genome_initialcondition_seeds.mat');
allgenomes = A;% allgenomes is a cell array containing all the seeds


load('Adj_asym_sig_seeds.mat'); % A2 is a column matrix of seeds for the random number generator for each value of p_sig ,p_asym and p_adj

% initial conditions tested here are all single steady states. not testing combinations.
numinit = 10;% number of initial conditions to be tried for each genome

prob1=0:0.1:1;% all possible values of p_sig, p_asym and p_adj

all_cellstates=dec2bin(0 : (2^N - 1),N)-'0';% all possible cell states with N node genomes
 


for genomeID = 1:length(allgenomes) 
      
    rand1 = str2double(allgenomes{genomeID});
  rand('state',rand1);% set seed for random number generator
  
  % generate the genome = basins + steady states
  B = randi(2^N,1,2^N);% basin
  uniqB = unique(B);
  numss = length(uniqB);
  
  genome1 = {cell(length(uniqB),1), cell(length(uniqB),1)};% genome = {{steady states}, {basins}}
  
  N_fp = 0;% number of fixed points
  N_os = 0;% number of oscillatory steady states
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
    
     if sum_ss == 1
         N_fp = N_fp + 1;
     else
         N_os = N_os + 1;
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
  
  
  fname = sprintf('lineage_maps_genome%d_N%d.mat', genomeID, N);
  load(fname);
    

    for initID = 1 : length(initC)
        count = 0;
        for Sig=1:length(prob1)
          for adj1=1:length(prob1)

            for asym1= 1:length(prob1)

		% write genome
		 fprintf(fid,'%d,',N);% number of genes

		  % write N_fp
		  fprintf(fid,'%d,',N_fp);

 		  % write N_os
 		  fprintf(fid,'%d,',N_os);

		  % write basin
		  for i1 = 1 : 2^N - 1
		      fprintf(fid,'%d;',B(i1));
		  end
		  fprintf(fid,'%d,',B(end));

		  %write ss
		   for i1 = 1 : 2^N - 1
		      fprintf(fid,'%d;',S(i1));
		   end
		  fprintf(fid,'%d,',S(end));

		   disp('written genome');




                  count = count + 1;

                  ruleID=A2(count);% seed for mersenne twister
                 
                  rand('state',(ruleID + rand1)/2);

                  sig=rand(1,N) <= prob1(Sig);% identities of signalling molecules
                  
                  % write sig
                  fprintf(fid,'%d,',prob1(Sig));
                  
                  % write sigvec
                  for i1 = 1 : length(sig) - 1
                      fprintf(fid,'%d;',sig(i1));
                  end
                  fprintf(fid,'%d,',sig(end));
                  
                  
		disp('written sig');

                  cell_adj = prob1(adj1);% probability that a given cell type receives signals from another given cell type
                  adjacencies = rand(2^N) <= cell_adj; % signaling adjacency matrix
                  adj2 = reshape(adjacencies',1,numel(adjacencies));
                  
                  % write adj
                  fprintf(fid,'%d,',cell_adj);
                  
                  % write adjmat
                  for i1 = 1 : length(adj2) - 1
                      fprintf(fid,'%d;',adj2(i1));
                  end
                  fprintf(fid,'%d,',adj2(end));
                  
                  disp('written adj');

                  p_asym=prob1(asym1); % probability of an 'on' gene being 'off' in one of the daughters
                  
                  % write asym
                  fprintf(fid,'%d,',p_asym);
                  
                  daughtermatrix = rand(2^N,N) >= p_asym; % p_asym = 0 implies daughter cells are identical to mother cells, p_asym = 1 implies daughter cells can only be [ 0 0 0 ...]. use this as a filter to get identities of daughter cells

                  mother_index = [];% identities of mother cells for each daughter cell produced 
                  Daughters = []; % set of daughters produced by all cell types present 
                  
                  divmat = zeros(2^N);% a binary matrix, where divmat(i,j) = 1 implies cell j is a daughter of cell i. only steady state rows are filled
                  
                  for i1 = 1 : numss
                    ss = genome1{1}{i1};
                    d1 = [];
                    for i2 = 1 : size(ss,1) % in the case of oscillations, all cell-types are considered to be present, therefore, daughters of all oscillating cells are produced at the same time. in this case, even for asym = 0, >1 daughter cells are produced
                      d1 = [d1;daughtermatrix * diag(ss(i2,:))];
                      
                      % find ss and d1. fill divmat
                      r = find(ismember(all_cellstates,ss(i2,:),'rows'));
                      c = ismember(all_cellstates,d1,'rows');
                      
                      divmat(r,c) = 1;
                    end
                  end
                  
                  divmat2 = reshape(divmat',1,numel(divmat));
                  
                  % write divmat
                  for i1 = 1 : length(divmat2) - 1
                      fprintf(fid,'%d;',divmat2(i1));
                  end
                  fprintf(fid,'%d,',divmat2(end));
 		
			disp('written asym');                 
                  
                  linmap1 = lineage_maps{initID,Sig,asym1,adj1};% edge list

                  if any(linmap1)
                 	 % write lineage map
                 	 linmap2 = reshape(linmap1',1,numel(linmap1));
                  
                  	 for i1 = 1 : length(linmap2) - 1
                     	 	fprintf(fid,'%d;',linmap2(i1));
                 	 end
                 	 fprintf(fid,'%d\n',linmap2(end));
                 	  
                    	 disp('written linmap');
		  else
			fprintf(fid,'%d\n',0);
			disp('no linmap here');
	          end
                                  
            end
          end

        end
    end
    
end

fclose(fid);

disp('closed file');
