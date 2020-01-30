% 13/8/2019

%write data into a text file. delimiters= ';' within a
% single object AND ',' between objects

N = 3;
fname = pwd;

for genomeID = 1:20
	filename = sprintf('trajectories_N3_genome%d.txt',genomeID);
	fid = fopen(filename,'a+');
	
	disp('opened file for reading and writing');
	Dat = textread(filename,'%s');% Dat is a cell array where each line is a separate cell, and it is read as a string
	disp('read file');
	dl1 = len(Dat);% number of trajectories already calculated
	T = cell(11,11,11,10);
	count = 0;
	for Sig = 1:11
		for asym1 = 1:11
			for adj1 = 1:11
				for initID = 1:10
					count = count + 1;
					% find trajectories for 0.1 of data
					if (rand < 0.1)&&(count < dl1)
						[traj,target] = explore_trajectories(N,genomeID,initID,Sig,asym1,adj1,fname);
						disp('trajectory code done')
						% write all IDs and parameters for each graph, number of nodes, 
						if isempty(traj{1})==0	
							for i1 = 1:length(traj)
								T1 = traj{i1};
								fprintf(fid,'%d,',N);% number of genes
								fprintf(fid,'%d,',genomeID);
								fprintf(fid,'%d,',initID);
								fprintf(fid,'%d,',Sig);
								fprintf(fid,'%d,',asym1);
								fprintf(fid,'%d,',adj1);
								fprintf(fid,'%d,',size(T1,1));% number of steps
								fprintf(fid,'%d,',target(i1));% whether the target graph is reached (regenerated or not)
								for i2 = 1:size(T1,1)
									fprintf(fid,'%3.2f,',T1(i2,1)/T1(i2,2));% fraction of cells in current state that are in target
								end
								for i2 = 1:size(T1,1)-1
									fprintf(fid,'%d,',T1(i2,2));% number of cells in current state
								end
								fprintf(fid,'%d\n',T1(end,2));
			
								disp('written results to file');
							end
						end
						% save all results in matfile
						%R1{initID,Sig,asym1,adj1} = r1;
						%R2{initID,Sig,asym1,adj1} = r2;
						%Frac_r1{initID,Sig,asym1,adj1} = frac_r1;
						%Frac_r2{initID,Sig,asym1,adj1} = frac_r2;
						%Numsteps1{initID,Sig,asym1,adj1} = numsteps1;
						%Numsteps2{initID,Sig,asym1,adj1} = numsteps2;
						%save(sprintf('reproducibility_data_genome%d_N%d.mat', genomeID, N),'R1','R2','Frac_r1','Frac_r2','Numsteps1','Numsteps2');
						%disp('saved results');
						disp([genomeID,initID,Sig,asym1,adj1]);
					end
 				end
			end
		end
	end
	fclose(fid);
	disp('closed file');
end

