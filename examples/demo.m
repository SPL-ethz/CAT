
clc;
PD = ProblemDefinition;
PD.sol_method = 'movingpivot';
tic
example_TemperatureProfile
f = toc;
fprintf('Solution Time for Moving Pivot Method %4.2f s\n',f)
PDMP = PD;
clear PD
PD = ProblemDefinition;
PD.sol_method = 'centraldifference';
tic
example_TemperatureProfile
f= toc;
fprintf('Solution Time for Central Differences Method %4.2f s\n',f)
PDCD = PD;
clear PD
PD = ProblemDefinition;
PD.sol_method = 'hires';
tic
example_TemperatureProfile
f=toc;
fprintf('Solution Time for High Resolution Method %4.2fs\n',f)
PDHR = PD;
clear PD

figure(1) % final PSD
subplot(1,2,1) % number weighted
plot(PDCD.calc_dist(end).y,PDCD.calc_dist(end).F,'k-o')
hold on
plot(PDMP.calc_dist(end).y,PDMP.calc_dist(end).F,'b-x')
plot(PDHR.calc_dist(end).y,PDHR.calc_dist(end).F,'r-^')
grid on
xlabel('L [\mum]')
ylabel('f [1/(g \mum)]')
legend('Central Differences','Moving Pivot','High Resolution')

subplot(1,2,2)
plot(PDCD.calc_dist(end).y,PDCD.calc_dist(end).F.*PDCD.calc_dist(end).y.^3,'k-o')
hold on
plot(PDMP.calc_dist(end).y,PDMP.calc_dist(end).F.*PDMP.calc_dist(end).y.^3,'b-x')
plot(PDHR.calc_dist(end).y,PDHR.calc_dist(end).F.*PDHR.calc_dist(end).y.^3,'r-^')
grid on
xlabel('L [\mum]')
ylabel('f L^3 [\mum^3/(g \mum)]')
legend('Central Differences','Moving Pivot','High Resolution')

figure(2) % number and volume weighted average sizes over time
subplot(1,2,1)
plot(PDCD.calc_time,moments(PDCD.calc_dist,1)./moments(PDCD.calc_dist,0),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,1)./moments(PDMP.calc_dist,0),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,1)./moments(PDHR.calc_dist,0),'r-^')
grid on
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu1/\mu0')

subplot(1,2,2)
plot(PDCD.calc_time,moments(PDCD.calc_dist,4)./moments(PDCD.calc_dist,3),'k-o')
hold on
plot(PDMP.calc_time,moments(PDMP.calc_dist,4)./moments(PDMP.calc_dist,3),'b-x')
plot(PDHR.calc_time,moments(PDHR.calc_dist,4)./moments(PDHR.calc_dist,3),'r-^')
axis tight
grid on
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
xlabel('Time [s]')
ylabel('\mu4/\mu3')

figure(3) % Phase diagram
Tcalc = PDMP.Tprofile(PDMP.calc_time);
Tvec = linspace(min(Tcalc)-5,max(Tcalc)+5);
plot(Tvec,PDMP.solubility(Tvec),'--','linewidth',1.5)
hold on
plot(PDCD.Tprofile(PDCD.calc_time),PDCD.calc_conc,'k-o')
plot(PDCD.Tprofile(PDMP.calc_time),PDMP.calc_conc,'b-x')
plot(PDCD.Tprofile(PDHR.calc_time),PDHR.calc_conc,'r-^')
axis tight
legend('Central Differences','Moving Pivot','High Resolution','location','southeast')
grid on
xlabel('Temperature [\circC]')
ylabel('Concentration [g/g]')