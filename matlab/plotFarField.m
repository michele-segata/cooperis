clear
clf

% incident direction the RIS should be configured to optimize
phiI = 0;
thetaI = 0;

% reflection direction the RIS should be configured to optimize
phiR = -45;
thetaR = 45;

% actual incident direction of the transmitter
thetaTX = 0;
phiTX = 0;

% frequency of operation, in Hz
f = 25e9;
% number of discrete phases that can be assigned to each element of the RIS
% for example, for n = 4, you have 4 possible phases: 0, pi/2, pi, 3pi/2
n = 8;
% number of elements per lambda. 5 means in one wavelength you have 5
% elements
cellsPerLambda = 5;
% size of one side of the RIS, in multiples of the wavelength
lambdaSize = 5;
% to compute the total number of elements (we assume a square surface):
% (cellsPerLambda * lambdaSize)^2

% angular resolution of the far field plot in degrees, for phi and theta
% by default it is 1 degree. Setting it to 0.1 degrees makes the plot
% smoother but increases the computation time
dp_deg = 1;
dt_deg = 1;

% compute the optimal phases given the desired incidence and reflection
% directions
[phases] = risCoding(phiR, thetaR, phiI, thetaI, f, n, cellsPerLambda, lambdaSize);
% compute the gain for all possible combinations of phi and theta (position
% of the receiver)
[power, phi, theta] = risFarField(phases, dp_deg, dt_deg, phiTX, thetaTX, f, cellsPerLambda, lambdaSize);

% plot the linear gain as a heatmap
h = surf(phi, theta, power);
% decomment the following line if you want to plot the phases assigned to
% the elements
% surf(phases)
fprintf("Maximum gain (linear): %f\n", max(max(power)));
fprintf("Maximum gain (dB)    : %f\n", 10*log10(max(max(power))));

% some graphical adjustments
set(h,'edgecolor','none')
xlabel("phi")
ylabel("theta")
colorbar;
xlim([-180 180]);
ylim([0 90]);
% uncomment this if you want to have a fixed color scale
% clim([0 105]);

% add an X where the maximum gain SHOULD BE located
% this might not coincide with the actual maximum, for example if the
% transmitter is placed in the wrong position
hold on;
ZMax=max(power(:));
plot3(phiR,thetaR,ZMax,'rx')

% view the heatmap from above
view([0 90])