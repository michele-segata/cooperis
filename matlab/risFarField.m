function [power, phi, theta] = risFarField(phases, dp_deg, dt_deg, phiTX, thetaTX, f, cellsPerLambda, lambdaSize)

[Fa, phi, theta] = risFarFieldPhases(phases, dp_deg, dt_deg, phiTX, thetaTX, f, cellsPerLambda, lambdaSize);

% computation of the far field
Dt = dt_deg * pi / 180;
Df = dp_deg * pi / 180;
A_theta = a_spherical(theta * pi / 180, Dt, Df);

% computing power
P = Fa.^2;

v2 = sum(P' * A_theta');
% computing normalized power in linear scale w.r.t. isotropic radiator
power = P * 2 * pi / v2;

end

