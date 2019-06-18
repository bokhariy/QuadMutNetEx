function [met]=QuadMutNetExMetrics_GENES_Struct(GenesNamesSolutionSet,mx,GenesNames,nets,C, k,np)
    Index = ismember(GenesNames, GenesNamesSolutionSet);
    sol =  find(Index==1);
    subNetworkSolution=nets(sol,sol);
    [met]= QuadMutNetExMetricsStruct(sol,mx,subNetworkSolution,C,k,np)