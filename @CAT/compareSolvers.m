function [CATS]=compareSolvers(kitty)

solverList={'movingpivot' 'centraldifference' 'hires'};

CATS(length(solverList)) = CAT;

for i = 1:length(CATS)
    kitloc = kitty;
    kitloc.sol_method = solverList{i};
    
    if ~strcmpi(solverList{i},kitty) || (strcmpi(solverList{i},kitty) && isempty(kitty.calc_time))
        kitloc.calc_time = [];
        kitloc.calc_dist = [];
        kitloc.calc_conc = [];
        kitloc.solve;
    end
    CATS(i) = kitloc;
end
    