function [F, phi, theta] = risFarFieldPhases(phases, dp_deg, dt_deg, phiTX, thetaTX, f, cellsPerLambda, lambdaSize)

c=299792458;
lambda = c/f;
k = 2 * pi / lambda;
du = lambda / cellsPerLambda;
kdu = k * du;
DN = lambdaSize * lambda;
DM = DN;
M = round(DM/du);
N = round(DN/du);

Ps = 360/dp_deg;
Ts = 90/dt_deg;
phi = linspace(-179,180,Ps);
theta = linspace(0,90,Ts+1);
phi_rad = phi * pi / 180;
theta_rad = theta * pi / 180;

F = 0;
for m=0:M-1
    for n=0:N-1

        % phase shift due to the position of the transmitter
        alpha = kdu * (n*sin(thetaTX)*cos(phiTX) + m*sin(thetaTX)*sin(phiTX));

        % overall phase shift considering: position of the transmitter
        % (alpha), coding (phases), and all possible positions of receivers
        phase = exp(-1j*( ...
            alpha + ...
            phases(m+1,n+1) - ...
            kdu * n.* sin(theta_rad)'* cos(phi_rad) - ...
            kdu * m.* sin(theta_rad)'* sin(phi_rad)) ...
        );

        F = F + phase;
    end
end

% after summing, for each theta and phi, all the phases of the signals
% coming from the different elements, we compute the absolute value (length
% of all the vectors), basically measuring constructive and distructive
% interference
F = abs(F);
end

