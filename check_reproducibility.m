% 2/11/2018
% reproducibility: at least one single cell-type of an organism should be capable of reproducing the complete organism.
% inputs: 
% 1. numss = number of steady states of the gene regulatory network
% 2. B: what basins do each of the different cell types belong to?
% 3. sig: identities of the molecules that act as signals
% 4. Daughters: the set of daughter cells produced by the different steady staes
% 5. mother_index: indicates which basin the mother cell of any given daughter cell belonged to 
% 6. adjacencies: a binary matrix indicating source and recipient cells for signals
% 7. s_final: set of steady cell types present in the final organism
% 8. N: number of genes in the system 

function r = check_reproducibility(numss,B,sig,Daughters,mother_index,adjacencies,s_final,N)
          
        r = 0;
        
        all_cellstates=dec2bin(0 : (2^N - 1),N)-'0';% all possible cell states with N node genomes
          
          
        for i1 = 1 : length(s_final)
          
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

            s1_id = s4_id;
            count1 = count1 + 1;
            
          end
          
          if isequal(s1_id,s_final)
            r = 1;
            break;
          end
        end
        
