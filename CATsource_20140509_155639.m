kitty = CAT; 

%%   Property: init_dist 
%   Initial particle size distribution 
%   Defined using Distribution class 
kitty.init_dist = Distribution([],[]);


%%   Property: init_seed 
%   Initial seed mass. 
%   Scalar value. The units must be consistent with those used for: 
%    * Initial concentration 
%    * Solubility function 
%    * Crystal density 
kitty.init_seed =[]; 


%%   Property: init_massmedium 
%   Initial mass of the continuous medium: total mass of solvent and 
%   antisolvent. 
%   Scalar value. The units must be consistent with those used for: 
%    * Initial concentration 
%    * Solubility function 
%    * Antisolvent profile 
%    * Nucleation rate 
kitty.init_massmedium =[]; 


%%   Property: init_conc 
%   Initial concentration. 
%   Scalar value, with units of mass solute per total solvent (continuous 
%   medium) mass. Use 'S=xx' for a solution with a defined Supersaturation. 
%   Units must be consistent with those used for: 
%    * Seed mass 
%    * Crystal density 
%    * Total initial solvent+antisolvent mass 
%    * Solubility function 
%    * Antisolvent profile 
%    * Nucleation rate 
kitty.init_conc =4; 


%%   Property: Tprofile 
%   Temperature profile: temperature as a function of time. 
%   Can be: anon. function, matrix or scalar. 
%   Units must be consistent with those used for: 
%    * Solubility function 
%    * Nucleation rate 
%    * Growth rate 
%    * Solubility 
kitty.Tprofile =[]; 


%%   Property: ASprofile 
%   Antisolvent profile: mass added as a function of time. 
%   Can be: anon. function, or matrix 
%   Units must be consistens with those used for: 
%    * Initial concentration 
%    * Solubility function 
%    * Nucleation rate 
%    * Growth rate 
%    * Solvent mass 
kitty.ASprofile =[]; 


%%   Property: solubility 
%   Solubility function as a function of temperature (T) and antisolvent mass fraction (xm). 
%   Can be an anon. function (Defined as: @(T,xm) or a scalar. 
%   Units must be mass solute per total solvent mass and consistent with: 
%    * Seed mass 
%    * Crystal density 
%    * Total initial solvent+antisolvent mass 
%    * Initial concentration 
%    * Antisolvent profile 
%    * Nucleation rate 
kitty.solubility =[]; 


%%   Property: rhoc 
%   Crystal density. 
%   Scalar value. The units must be consistent with those used for: 
%    * Seed mass 
%    * Initial concentration 
%    * Solubility function 
kitty.rhoc =x; 


%%   Property 
%   Shape factor 
%   Dimensionless scalar value, between 0 and 1. 
kitty.kv =[]; 


%%   Property: growthrate 
%   Growth rate as a function of supersaturation (S), temperature (T) and size (y). 
%   Defined as an anon. function: @(S,T,y) 
%   Units must be consistent with those used for: 
%    * Initial distribution 
%    * Temperature profile 
%    * Antisolvent profile 
%    * Nucleation rate 
%    * Solution time 
%   
%   The function should return a vector the same size as y 
kitty.growthrate =[]; 


%%   Property: nucleationrate 
%   Nucleation rate function 
%   Nucleation rate as a function of supersaturation (S), temperature (T) or (moments of) the distribution F. 
%   Defined as an anon. function: @(S,T,F) 
%   Units must be consistent with those used for: 
%    * Initial distribution 
%    * Temperature profile 
%    * Antisolvent profile 
%    * Growth rate 
%    * Solution time 
%   
%   The function should return a scalar value 
kitty.nucleationrate =[]; 


%%   Property: sol_time 
%   Solution time 
%   Vector or scalar. 
%   The units must be consistent with those used for: 
%    * Growth rate 
%    * Nucleation rate 
%    * Antisolvent profile 
%    * Temperature profile 
kitty.sol_time =[]; 


%%   Property: sol_method 
%   Solution method. 
%   Defines which numerical method to use. 
%   Use the solutionMethods function to see a list of available 
%   solvers. 
kitty.sol_method ='centraldifference'; 


%%   Property: sol_options 
%   Solver options 
kitty.sol_options =[]; 


