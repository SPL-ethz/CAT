function [f] = hiResNucleation(c0,cs0,f0,Dt,finput,Dx,xp,ndim,xp1_arr,xp2_arr,xp3_arr)

%% Addition of Nucleation term

mu3         =   sum(f0(3:end-1) .*Dx.dim1 .* xp1_arr.^3);

B         =   (finput.data.kinpar.kb*mu3*((c0-cs0)/cs0)^finput.data.kinpar.b)*Dt;
f = f0;
f(3)    =   f0(3)        +   B;
