clear;
clf;

phiI = 0;
thetaI = 0;
thetaTX = thetaI;
phiTX = phiI;
thetaR = 45;

f = 25e9;
n = 8;
cellsPerLambda = 5;
lambdaSize = 5;

dp_deg = 0.1;
dt_deg = 0.1;

phiRs = [-45 45 135];
thetaRs = [45 45 45];
average = false;
random = true;
seedn = 1;
autorange = true;

plotBeamSplit(phiRs, thetaRs, phiI, thetaI, phiTX, thetaTX, dp_deg, dt_deg, f, n, cellsPerLambda, lambdaSize, average, random, seedn, autorange, "test");
