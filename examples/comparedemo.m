
function cumparedemo(example)

global demo

if nargin == 0
%    example = 'example_ASandTProfile'; 
   example = 'example_Nucleation'; 
   nBinsv = [200 500 500]; %number of bins for MP,CD,HR respectively;
end

close all;clc;
addALLthepaths
PD = ProblemDefinition;
PD.sol_method = 'movingpivot';
nBins = nBinsv(1);
tic
eval(example)
f = toc;
fprintf('Solution Time for Moving Pivot Method %4.2f s\n',f)
PDMP = PD;
clear PD
PD = ProblemDefinition;
PD.sol_method = 'centraldifference';
nBins = nBinsv(2);
tic
eval(example)
f= toc;
fprintf('Solution Time for Central Differences Method %4.2f s\n',f)
PDCD = PD;
clear PD
PD = ProblemDefinition;
PD.sol_method = 'hires';
nBins = nBinsv(3);
tic
eval(example)
f=toc;
fprintf('Solution Time for High Resolution Method %4.2fs\n',f)
PDHR = PD;
clear PD

figure(1) % final PSD
subplot(1,2,1) % number weighted
plot(PDCD.calc_dist(1).y,PDCD.calc_dist(1).F,'g--')
hold on
plot(PDCD.calc_dist(end).y,PDCD.calc_dist(end).F,'k-o')
plot(PDMP.calc_dist(end).y,PDMP.calc_dist(end).F,'b-x')
plot(PDHR.calc_dist(end).y,PDHR.calc_dist(end).F,'r-^')
grid on
xlabel('L [\mum]')
ylabel('f [1/(g \mum)]')
legend('Initial','Central Differences','Moving Pivot','High Resolution')

subplot(1,2,2)
plot(PDCD.calc_dist(1).y,PDCD.calc_dist(1).F.*PDCD.calc_dist(1).y.^3,'g--')
hold on
plot(PDCD.calc_dist(end).y,PDCD.calc_dist(end).F.*PDCD.calc_dist(end).y.^3,'k-o')
plot(PDMP.calc_dist(end).y,PDMP.calc_dist(end).F.*PDMP.calc_dist(end).y.^3,'b-x')
plot(PDHR.calc_dist(end).y,PDHR.calc_dist(end).F.*PDHR.calc_dist(end).y.^3,'r-^')
grid on
xlabel('L [\mum]')
ylabel('f L^3 [\mum^3/(g \mum)]')
legend('Initial','Central Differences','Moving Pivot','High Resolution')

figure(2) % number and volume weighted average sizes over time
subplot(2,2,1)
plot(PDCD.calc_time,moments(PDCD.calc_dist,1)./moments(PDCD.calc_dist,0),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,1)./moments(PDMP.calc_dist,0),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,1)./moments(PDHR.calc_dist,0),'r-^')
grid on
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu1/\mu0')

subplot(2,2,2)
plot(PDCD.calc_time,moments(PDCD.calc_dist,4)./moments(PDCD.calc_dist,3),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,4)./moments(PDMP.calc_dist,3),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,4)./moments(PDHR.calc_dist,3),'r-^')
axis tight
grid on
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu4/\mu3')

subplot(2,2,3)
plot(PDCD.calc_time,moments(PDCD.calc_dist,0),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,0),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,0),'r-^')
axis tight
grid on
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu0')

subplot(2,2,4)
plot(PDCD.calc_time,moments(PDCD.calc_dist,3),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,3),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,3),'r-^')
axis tight
grid on
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu3')

figure(3) % Process
subplot(2,2,1)
plot(PDCD.calc_time,PDCD.calc_conc,'k-o','linewidth',1.5)
hold on
plot(PDMP.calc_time,PDMP.calc_conc,'b-x','linewidth',1.5)
plot(PDHR.calc_time,PDHR.calc_conc,'r-^','linewidth',1.5)
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Time [s]')
ylabel('Concentration [g/g]')

subplot(2,2,2)
plot(PDCD.calc_time,PDCD.calc_conc(:)./PDCD.solubility(PDCD.Tprofile(PDCD.calc_time(:))),'k-o','linewidth',1.5)
hold on
plot(PDMP.calc_time,PDMP.calc_conc./PDMP.solubility(PDMP.Tprofile(PDMP.calc_time)),'b-x','linewidth',1.5)
plot(PDHR.calc_time,PDHR.calc_conc./PDHR.solubility(PDHR.Tprofile(PDHR.calc_time)),'r-^','linewidth',1.5)
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Time [s]')
ylabel('Supersaturation [-]')

subplot(2,2,3)
plot(PDCD.calc_time,PDCD.Tprofile(PDCD.calc_time),'k-o','linewidth',1.5)
hold on
plot(PDMP.calc_time,PDMP.Tprofile(PDMP.calc_time),'b-x','linewidth',1.5)
plot(PDHR.calc_time,PDHR.Tprofile(PDHR.calc_time),'r-^','linewidth',1.5)
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Time [s]')
ylabel('Temperature [\circC]')

subplot(2,2,4)
mscalc = PDCD.init_massmedium + PDCD.ASprofile(PDCD.calc_time)-PDCD.ASprofile(0);
plot(PDCD.calc_time,mscalc,'k-o','linewidth',1.5)
hold on
mscalc = PDMP.init_massmedium + PDMP.ASprofile(PDMP.calc_time)-PDMP.ASprofile(0);
plot(PDMP.calc_time,mscalc,'b-x','linewidth',1.5)
mscalc = PDHR.init_massmedium + PDHR.ASprofile(PDHR.calc_time)-PDHR.ASprofile(0);
plot(PDHR.calc_time,mscalc,'r-^','linewidth',1.5)
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Time [s]')
ylabel('Total mass Solvent + Antisolvent [g]')

figure(4) % Phase diagram
Tvec = linspace(min(PDMP.Tprofile(PDMP.calc_time))-5,max(PDMP.Tprofile(PDMP.calc_time))+5);
plot(Tvec,PDMP.solubility(Tvec),'--','linewidth',1.5)
hold on
plot(PDCD.Tprofile(PDCD.calc_time),PDCD.calc_conc,'k-o')
plot(PDCD.Tprofile(PDMP.calc_time),PDMP.calc_conc,'b-x')
plot(PDCD.Tprofile(PDHR.calc_time),PDHR.calc_conc,'r-^')
axis tight
legend('','Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Temperature [\circC]')
ylabel('Concentration [g/g]')