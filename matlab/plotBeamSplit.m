function [] = plotBeamSplit(phiRs, thetaRs, phiI, thetaI, phiTX, thetaTX, dp_deg, dt_deg, f, n, cellsPerLambda, lambdaSize, average, random, seedn, autorange, filename)

if average && random
    return
end

titlePhis = "phiRs = [";
titleThetas = "thetaRs = [";
for i = 1:size(phiRs, 2)
    [phases] = risCoding(phiRs(i), thetaRs(i), phiI, thetaI, f, n, cellsPerLambda, lambdaSize);
    PHAs{i} = phases;
    AMPs{i} = ones(size(phases));
    titlePhis = strcat(titlePhis, num2str(phiRs(i)));
    titleThetas = strcat(titleThetas, num2str(thetaRs(i)));
    if i ~= size(phiRs, 2)
        titlePhis = strcat(titlePhis, ", ");
        titleThetas = strcat(titleThetas, ", ");
    end
end
params = strcat(titleThetas, "], ", titlePhis, "], strategy=");
if average
    params = strcat(params, "average");
else
    if random
        params = strcat(params, "random");
    else
        params = strcat(params, "mode");
    end
end

step = 2*pi/n;
[~, pha] = CFG_mode_per_cell(AMPs, PHAs, 2*pi-step, step, 0.1, average, random, seedn);

if (average)
    par_phases_file = strcat( ...
        "phiR_", regexprep(num2str(phiRs), "  *", "_"), "_thetaR_", regexprep(num2str(thetaRs), "  *", "_"), ...
        "_phiI_", num2str(phiI), "_thetaI_", num2str(thetaI), ...
        "_phiTX_", num2str(phiTX), "_thetaTX_", num2str(thetaTX), ...
        "_n_", num2str(n), "_pl_", num2str(cellsPerLambda), "_nl_", num2str(lambdaSize));
    writematrix(pha, strcat("phases_", par_phases_file, ".csv"));
end

[power, phi, theta] = risFarField(pha, dp_deg, dt_deg, phiTX, thetaTX, f, cellsPerLambda, lambdaSize);

h = surf(phi, theta, power);
ax = gca;
ax.FontSize = 10;
title(params);
max(max(power))
set(h,'edgecolor','none')
xlabel("phi")
ylabel("theta")
colorbar;
if ~autorange
    clim([0 110]);
end
xlim([-180 180]);
ylim([0 90]);

view([0 90])

saveas(gcf, filename);

end