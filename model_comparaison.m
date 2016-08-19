clear all

% Simulate the data

model_experiment.type = 'first_order_homogenous'; % first_order_homogenous, first_order, second order, Navier-Stokes
model_experiment.coefficients = [1, 2]; % 1, 2, 3, ... coefficients
coefficients=[1, 0.1, 0, 1, 2]; %  [y_initial, delta_t, t(0), t(end)]
noise.isnoise = true;
noise.sigma_noise = 0.05;


q_experiment=makedist('Normal','mu',0,'sigma',noise.sigma_noise); % distribution normale centrée en 0

measurements = simulate_data(model_experiment, coefficients, noise, q_experiment); % measurements.y, measurements.t, measurements.number

plot(measurements.t,measurements.y)
%% Computation of the model


model = 'first_order_homogenous';
prior.lower_a = 0; % borne inférieure pour la prior uniforme sur a
prior.upper_a = 2; % borne supérieure
prior.lower_b = 0; % idem pour b
prior.upper_b = 4;


number_samples = 100; % nombres de d'échantillons
number_iteration = 100; % nombre d'itérations pour la chaîne de Markov
sigma_algorithm=0.1;

scheme = randi([-1 1],1,number_samples); % upwind 1, downwind -1, center 0
several_scheme = true;

[ parameters ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);
[ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, several_scheme, number_iteration );

figure()
histogram(parameters.scheme)
title('Distribution of the schemes');
xlabel('Downwind, Center, Upwind');


if strcmp(model,'first_order_homogenous')
    a_upwind=parameters.coefficients_upwind(1,:);
    a_downwind=parameters.coefficients_downwind(1,:);
    a_center=parameters.coefficients_center(1,:);
    
    figure()
    histogram(a_upwind,'BinWidth', noise.sigma_noise)
    title('Upwind');
    xlabel('a');
    
    figure()
    histogram(a_center,'BinWidth', noise.sigma_noise)
    title('Center');
    xlabel('a');
    
    figure()
    histogram(a_downwind,'BinWidth', noise.sigma_noise)
    title('Downwind');
    xlabel('a');
elseif strcmp(model,'first_order')
    
    a_upwind=parameters.coefficients_upwind(1,:);
    a_downwind=parameters.coefficients_downwind(1,:);
    a_center=parameters.coefficients_center(1,:);
    b_upwind=parameters.coefficients_upwind(2,:);
    b_downwind=parameters.coefficients_downwind(2,:);
    b_center=parameters.coefficients_center(2,:);    
    
    figure()
    histogram(a_upwind,'BinWidth', noise.sigma_noise)
    title(parameters.scheme);
    xlabel('a');    
    figure()
    histogram(b_upwind,'BinWidth', noise.sigma_noise)
    title('Upwind');
    xlabel('b');   
    
    figure()
    histogram(a_center,'BinWidth', noise.sigma_noise)
    title('Center');
    xlabel('a');    
    figure()
    histogram(b_center,'BinWidth', noise.sigma_noise)
    title('Center');
    xlabel('b');
    
    figure()
    histogram(a_downwind,'BinWidth', noise.sigma_noise)
    title('Downwind');
    xlabel('a');    
    figure()
    histogram(b_downwind,'BinWidth', noise.sigma_noise)
    title('Downwind');
    xlabel('b');    
end


    
    
    