QuadMutNetEx: Quadratic Driver Mutation Networks Explorer
Version: 1.0
Date: June-17-2019
contact: bokhariya@vcu.edu tarodz@vcu.edu

---------------------------------------------------------------------------------------------------------------------------------------------
Requirements:
1. Matlab
2- GUROBI OPTIMIZATION Version:7
Note: Matlab should be directed to Gurobi directory
----------------------------------------------------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------------------------------------------------
USAGE:

*utlRoulette, selectSubsetVer and Gurobi_optimization_nets are function needed and used in QuadMutNetEx

A- QuadMutNetEx:

	INPUT :
		Describtion of QuadMutNetEx arguments:
		1- G: 
		 	Mutation matrix in a matlab sparse format. Patients have a serial number from 1 to m. Genes have serial number from 1 to n.

		2- GenesNames:
			vector of genes names orderd as in the G sparse matrix(from 1 to n).

        3- nets:
            a matrix of one or multiple binary networks added together.

		4- t:
			Integer number of iterations desired.

		5- n:
		 	Integer number of genes to be included in the submatrix G' each iteration. Each iteration a matrix G' is chosen rendomly 
		 	from the big matrix G.

		6- c:
		 	Parameter c (As c gets bigger as the number of genes in solution shrinks).

		7- k:
			Parameter k (High value of k panalize for overlap).

        8- np:
            Parameter np (High value of np increases the effect of the network(s)

		9- GenesSelProb: (optional)
			vector of probabilities of each gene to be chosen randomly.The probabilities should be in one column orderd as 
			created serial number in the G_file(from 1 to n). The sum of the 
			colmn should be equal 1.

			Note: If left empty then each gene will have equal chance to be chosen
		

	OUTPUT:
		
		objective_function: Solution score that is the objective function we want to minimize.
		x: Selected genes numbers as given and orderd in the sparse matrix
		solutionSet: names of selected genes
        subNetSol: SolutionSet submatrix of summations of matrix(es) nets.


see also QuadMutNetEx_Example.m