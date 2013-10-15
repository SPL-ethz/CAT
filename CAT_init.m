function CAT_init

disp('Initializing...')

startstr = pwd;
addpath(startstr)
addpath(strcat(startstr,'/general'))
addpath(strcat(startstr,'/solvers'))
addpath(strcat(startstr,'/solvers/hires'))
addpath(strcat(startstr,'/solvers/movingpivot'))
addpath(strcat(startstr,'/solvers/centraldifference'))

disp('Done.')