function [objective_function, x, solutionSet, subNetSol] = QuadMutNetEx(G, GenesNames, nets, t, n, c, k, np, GenesSelProb)
%
%Yahya Bokhari  June-13-2019
%
%
% matlab dirctory should be where Gurobi sofrware saved (ex: /opt/gurobi702/linux64/matlab/)
%
%
% Describtion of QuadMutNetEx_FromFiles arguments:
% 1- G: 
% 	mutation matrix in a sparse format. Patients have a serial number from 1 to m. Genes have serial number from 1 to n.
% 2- GenesNames:
%	vector of genes orderd as in the G sparse matrix(from 1 to n).
% 3- nets:
%	a matrix of one or multiple binary networks added together.
% 4- t:
%	Number of iterations desired.
% 5- n:
% 	Number of genes to be included in the submatrix G' each iteration.
% 6- c:
% 	Parameter c (As c gets bigger as the number of genes in solution shrinks).
% 7- k:
%	Parameter k (High value of k panalize for overlap).
% 8- np:
%   Parameter np (High value of np increases the effect of the network(s)
%   in the solution set.
% 9- GenesSelProb: (optional)
%	vector of probabilities of each gene to be chosen randomly. The probabilities should be in one column orderd as 
%	created serial number in the G_file(from 1 to n). The sum of the colmn should be equal 1.
%	
if (nargin<8) 
 		disp('Usage: QuadMutNetEx(G,GenesNames,t,n,c,k,GenesSelProb)')
 		disp(' G: Mutatation sparse matrix.') 
 		disp(' GenesNames: a vector contains list of genes orderd as in the G sparse matrix..')
        disp(' nets: a matrix of one or multiple binary networks added together.')
 		disp(' t: Number of iterations.')
 		disp(' n: Number of genes to be included in the submatrix G''')
 		disp(' c: Parameter c.')
 		disp(' k: Parameter k.')
        disp(' np: Parameter np.')
 		disp(' GenesSelProb: a vector of Genes selection probabilities.')
 		disp(' (Type: help QuadMutNetEx OR see README.txt for more details).')
 	return;
 end 

mySparseMatrix=G;
[sCnt,gCnt]= size(mySparseMatrix);

if (nargin<9)
	selection_prob=[1:length(GenesNames)]./length(GenesNames);
else
	selection_prob=GenesSelProb;
end


netsSumF = np * nets;

qObjPrev=10000000;
bestQ=10000000;  % we want to minimize
best_sol=[];

prev_x = zeros(1,n);  % n is the number of genes to be tried each iteration.
prev_genes = [];

Temp=2;

obs=zeros(1,t);
obsB=zeros(1,t);
tic

for j = 1 :t
    
	[curr_genes] =selectSubsetVer(prev_x, prev_genes, selection_prob);

	curr_G=mySparseMatrix(:,curr_genes);  %  new gene names and G_submatrix
	subNet=netsSumF(curr_genes,curr_genes);
    

	[curr_x,qObjNoN,R]= Gurobi_optimization_nets(curr_G,subNet,c,k);  %Gurobi part
    	curr_sol=curr_genes(find(curr_x>0.5));
    

	qObj= sCnt+qObjNoN;
    
    
	if (qObj<bestQ)
        %store to be returned
		bestQ=qObj;
        best_sol=curr_sol;
       
        
	end
	obs(j)=qObj;
	obsB(j)=bestQ;

	deltaQ =  qObj  - qObjPrev;  
	luck = exp(-(deltaQ/Temp));
	coin= rand(1);

	if (coin <= luck)	%accept
		qObjPrev=qObj;
        prev_genes=curr_genes;
        prev_x=curr_x;        
	end  

	
end  % end of iter

%%############  RESULTS  ###############


objective_function=bestQ;
x=best_sol;
solutionSet=GenesNames(x);
subNetSol=nets(best_sol,best_sol); % to be used in mitrics
toc;