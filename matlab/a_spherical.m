% computes the size of a spherical element for a given elevation theta and
% the angular resolutions for theta and phi
function [a] = a_spherical(theta, dt, dp)
    a = zeros(size(theta));
    a(theta == 0) = dp * (1 - cos(dt/2));
    a(theta == pi/2) = dp * (cos(theta(theta == pi/2)-dt/2));
    a(theta > 0 & theta < pi/2) = dp * (cos(theta(theta > 0 & theta < pi/2)-dt/2) - cos(theta(theta > 0 & theta < pi/2)+dt/2));
end