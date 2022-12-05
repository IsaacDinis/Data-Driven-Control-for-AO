%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHILIPPE SCHUCHERT            %
% SCI-STI-AK, EPFL              %
% philippe.schuchert@epfl.ch    %
% March 2021                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Iteratively solve the datadriven problem
%

function [controller, sol] =  datadriven(system,obj,cons,params,solver)
fprintf("\n  iter  |   slack    |     obj     |  decrease   | Total SOLVE time\n")

if nargin == 5
    params.solver = solver;
else
    params.solver = '';
end

switch params.solver
    case 'fusion'
        fprintf(" --------------- SOLVING USING MOSEK + MOSEK FUSION ----------------\n")
    otherwise
        fprintf(" ---------------------- SOLVING USING YALMIP ------------------------\n")
end
iter = 0; op = NaN; sol = [];
solveTime = 0;
while iter < params.maxIter
    iter = iter + 1;
    switch params.solver
        case 'fusion'
            
            [system.controller,sol,d] = utils.solveddMOSEK(system,obj,cons,params,sol);
        otherwise
            [system.controller,sol,d] = utils.solveddYALMIP(system,obj,cons,params,sol);
    end
    solveTime = solveTime + d.solvertime; % Does not include yalmiptime!
    if isnan(op)
        diffString = '           ';
        
    else
        diffString = num2str(sol.obj-op, '%.04e');
        
    end
    if sol.satisfyConstraints
        solString = num2str(abs(sol.obj),'%+.04e');
    else
        solString = '           ';
    end
    
    
    fprintf("   %03d  | %.04e | %s | %s | %8.04e \n",sol.nIter, abs(sol.slack), solString , diffString,solveTime)
    
    if sol.obj-op >= -params.tol
        break
    end
    if sol.satisfyConstraints
        op = sol.obj;
    end
    if sol.slack < 1e-4
        sol.satisfyConstraints  = 1;
    end
end

fprintf(" --------------------------- END SOLVE -----------------------------\n\n")

controller = system.controller;
sol = rmfield(sol,'satisfyConstraints');
end





