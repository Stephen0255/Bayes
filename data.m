% Simulate the data

model_experiment='first_order'; % first_order_homogenous, first_order, second order, Navier-Stokes
coefficients=[1, 0.1, 0, 1, 2]; %  [y_initial, delta_t, t(0), t(end)]
noise.isnoise = false;
noise.sigma_noise = 0.01;


q=makedist('Normal','mu',0,'sigma',noise.sigma_noise); % distribution normale centrée en 0

function [ measurements ] = simulate_data( model, coefficients, noise, q )
% Generate the experimental data
%   Choose the model : first order equation, second order, ...
%   Choose the coefficients of the equation and the level of noise (experimental accuracy)
%   The outpout will be the position.

measurements.y_initial = coefficients(1);
measurements.delta = coefficients(2);
t = coefficients(3):measurements.delta:coefficients(4);

if strcmp(model, 'first_order_homogenous')
    y = measurements.y_initial*exp(-a*t);
elseif strcmp(model, 'first_order')
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

measurements = simulate_data(model_experiment, coefficients, noise, q); % measurements.y, measurements.t, measurements.number

