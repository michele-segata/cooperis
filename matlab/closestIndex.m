function [index]= closestIndex(vector, value)
    % find the available value in the set of the available ones that is
    % closer to the target
    distance = mod(vector - value + pi, 2*pi) - pi;
    itofix = find(distance < -pi);
    distance(itofix) = distance(itofix) + 2*pi;
    distance = abs(distance);
    [~, index] = min(distance);
end