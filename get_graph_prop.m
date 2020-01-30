% 15/10/2018
% classify maximal lineage maps
% cyclic, chains, tree, other DAG

% separate all maps into connected components

%{
(remove self edges)
chains:
there is exactly one source node and one target node.
all other nodes have exactly one incoming and one outgoing edge
chain length= number of nodes

divergent trees (that are not chains): 
There is exactly one source node, but multiple target nodes 
all nodes, except the source node, receive exactly one incoming edge

convergent trees:
inverting all edges gives a divergent tree

other DAGs: none of the above three categories
%}

% We have used graph analysis codes from MIT strategic engine's graph_theory_toolbox, which is published under the BSD liscence (available at: http://strategic.mit.edu/downloads.php?page=matlab_networks)



function Tab1 = get_graph_prop(g1_list)

addpath('graph_theory_toolbox_MIT');% You will need to obtain the graph_theory toolbox from the source (http://strategic.mit.edu/downloads.php?page=matlab_networks)


if any(g1_list)
        Tab1 = [];% 11 col = [numcomponents, componentID, numnodes, numedges, number of source nodes, number of target nodes, cycle, chain, divergent tree, convergent tree, other DAG]
     
        g1 = zeros(max(max(g1_list)));
        g1(sub2ind(size(g1), g1_list(:,1), g1_list(:,2))) = 1;
        
        no_nodes = (sum(g1,1) + sum(g1,2)') == 0; % nodes with no input or output edges. remove these
        g1(no_nodes,:) = [];
        g1(:,no_nodes) = [];

        g2 = g1-diag(diag(g1));% remove self-edges
        g2 = (g2+g2'+eye(length(g2)))>0;% make undirected
        g2 = findConnComp(g2);% each row is a connected component
        Tab1(1) = length(g2);
        Tab1 = repmat(Tab1, size(g2,1), 1);

        for i4 = 1 : length(g2)
          Tab1(i4,2) = i4;
          
          c1 = g1(g2{i4},:);
          c1 = c1(:,g2{i4});
          c11=c1-diag(diag(c1));
          
          scc_c11 = tarjan(adj2adjL(c1));
          cyclic = length(scc_c11) < length(c11);
      
          % classify c1
          s1=sum(sum(c11)==0);% number of source nodes
          s2=sum(sum(c11')==0);% number of target nodes
          s3=max(sum(c11));% max in degree
          s4=max(sum(c11'));% max out degree
          Tab1(i4,[3,4,5,6]) = [length(c1), sum(sum(c1)), s1, s2];
          Tab1(i4,7) = double((s1 == 0)&&(s2 == 0));% cycle (single strongly connected component)
          Tab1(i4,8)=double((s1==1)&&(s2==1)&&(s3==1)&&(s4==1))*length(c1);% chain
          Tab1(i4,9)=double((s1==1)&&(s2>1)&&(s3==1)&&(s4>1))*length(c1);%divergent tree
          Tab1(i4,10)=double((s1>1)&&(s2==1)&&(s3>1)&(s4==1))*length(c1);%convergent tree
          Tab1(i4,11) = isequal(Tab1(i4,[7,8,9,10]),[0,0,0,0]) & (cyclic == 0);

        end
else
        Tab1 = zeros(1,11);
end



