function [ measurements ] = simulate_data( model, coefficients, noise, q )
% Generate the experimental data
%   Choose the model : first order equation, second order, ...
%   Choose the coefficients of the equation and the level of noise (experimental accuracy)
%   The outpout will be the position.

measurements.y_initial = coefficients(1);
measurements.delta = coefficients(2);
t = coefficients(3):coefficients(2):coefficients(4);

if strcmp(model.type, 'first_order_homogenous')
    a = model.coefficients(1);
    y = measurements.y_initial*exp(-a*t);
elseif strcmp(model.type, 'first_order')
    a = model.coefficients(1);
    b = model.coefficients(2);
    y = (measurements.y_initial+b/a)*exp(-a*t)-b/a;
end

measurements.y = y;
measurements.t = t;
measurements.number = length(measurements.y);

if noise.isnoise
    measurements.y = measurements.y + random(q);
    measurements.t = measurements.t + random(q);
end

end

