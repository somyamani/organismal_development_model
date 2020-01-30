% 15/6/2019

%write data into a text file. delimiters= ';' within a
% single object AND ',' between objects

N = 3;% number of genes
filename = 'reproducibility_N3_1.txt';
fid = fopen(filename,'a');

disp('opened file for writing');
fname = pwd;

for genomeID = 1:20
	
	R1 = cell(10,10,10,10);
	R2 = cell(10,10,10,10);
	Frac_r1 = cell(10,10,10,10);
	Frac_r2 = cell(10,10,10,10);
	Numsteps1 = cell(10,10,10,10);
	Numsteps2 = cell(10,10,10,10);
	
	for Sig = 1:10
		for asym1 = 1:10
			for adj1 = 1:10
				for initID = 1:10

					[r1,frac_r1,numsteps1,r2,frac_r2,numsteps2] = explore_reproducibility(N,genomeID,initID,Sig,asym1,adj1,fname);
					disp('repoducibility code done')
					% write all IDs and parameters for each graph, number of nodes, frac_r1, frac_r2
					if r1 ~= 2
						fprintf(fid,'%d,',N);% number of genes
						fprintf(fid,'%d,',genomeID);
						fprintf(fid,'%d,',initID);
						fprintf(fid,'%d,',Sig);
						fprintf(fid,'%d,',asym1);
						fprintf(fid,'%d,',adj1);
						fprintf(fid,'%d,',length(r1));% number of nodes in graph
						fprintf(fid,'%d,',frac_r1);% fraction of 'pluripotent' cells
						if r2 ~= 2
							fprintf(fid,'%d,',length(r2));% number of non-graph cells
						else
							fprintf(fid,'0,');% number of non-graph cells
						end
						fprintf(fid,'%d\n',frac_r2);% fraction of non-graph cells that map to input graph
						disp('written results to file');
					end
					% save all results in matfile
					R1{initID,Sig,asym1,adj1} = r1;
					R2{initID,Sig,asym1,adj1} = r2;
					Frac_r1{initID,Sig,asym1,adj1} = frac_r1;
					Frac_r2{initID,Sig,asym1,adj1} = frac_r2;
					Numsteps1{initID,Sig,asym1,adj1} = numsteps1;
					Numsteps2{initID,Sig,asym1,adj1} = numsteps2;
					save(sprintf('reproducibility_data_genome%d_N%d.mat', genomeID, N),'R1','R2','Frac_r1','Frac_r2','Numsteps1','Numsteps2');
					disp('saved results');
					disp([genomeID,initID,Sig,asym1,adj1]);
 				end
			end
		end
	end
end

fclose(fid);
disp('closed file');
