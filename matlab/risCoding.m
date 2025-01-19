function [phases, M, N] = risCoding(phiR, thetaR, phiI, thetaI, f, n, cellsPerLambda, lambdaSize)

c=299792458;
lambda = c/f;
k = 2 * pi / lambda;
du = lambda / cellsPerLambda;
kdu = k * du;
DN = lambdaSize * lambda;
DM = DN;
M = round(DM/du);
N = round(DN/du);

E = 2*pi / n * linspace(0, n-1, n);

phiR = phiR * pi / 180;
thetaR = thetaR * pi / 180;
thetaI = thetaI * pi / 180;
phiI = phiI * pi / 180;



% -sin(thetaI) * cos(phiI) removes the phase shift due to arrival on x
% -sin(thetaI) * sin(phiI) removes the phase shift due to arrival on y
% cos(phiR) * sin(thetaR) adds the phase shift for the reflection on x
% sin(phiR) * sin(thetaR) adds the phase shift for the reflection on y
dsx = kdu * (cos(phiR)*sin(thetaR)-sin(thetaI)*cos(phiI));%%
dsy = kdu * (sin(phiR)*sin(thetaR)-sin(thetaI)*sin(phiI));%%

PHI = zeros(M, N); % coding matrix
states = zeros(M, N);

for m=0:M-1
    for n=0:N-1
        PHI(m+1,n+1) = ((n*dsx) + (m*dsy));

        PHI(m+1,n+1) = mod(PHI(m+1,n+1), 2*pi);

        % find the available phase in the set of the available ones that is
        % closer to the target one
        p = closestIndex(E, PHI(m+1, n+1));
        states(m+1,n+1) = p;
        PHI(m+1,n+1) = E(p);
    end
end

phases = PHI;

end

