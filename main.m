clear all

% Simulate the data

model_experiment.type = 'first_order_homogenous'; % first_order_homogenous, first_order, second order, Navier-Stokes
model_experiment.coefficients = [1, 2]; % 1, 2, 3, ... coefficients
coefficients=[1, 0.1, 0, 1, 2]; %  [y_initial, delta_t, t(0), t(end)]
noise.isnoise = true;
noise.sigma = 0.05;

q_experiment=makedist('Normal','mu',0,'sigma',noise.sigma); % distribution normale centrée en 0

measurements = simulate_data(model_experiment, coefficients, noise, q_experiment); % measurements.y, measurements.t, measurements.number

plot(measurements.t,measurements.y)
%% Computation of the model

model = 'first_order_homogenous';

number_samples = 50; % nombres de d'échantillons
number_iteration = 100; % nombre d'itérations pour la chaîne de Markov
sigma_algorithm=0.1;

scheme_comparaison = true; 

if scheme_comparaison
    scheme = randi([-1 1],1,number_samples); % upwind 1, downwind -1, center 0
    several_scheme = true;
else
    scheme = 1; % upwind 1, downwind -1, center 0
    several_scheme = false;
end

prior.lower_a = 0; % borne inférieure pour la prior uniforme sur a
prior.upper_a = 2; % borne supérieure
prior.lower_b = 0; % idem pour b
prior.upper_b = 4;

[ parameters ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);

[ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, several_scheme, number_iteration );

create_histogram( model, parameters, several_scheme, noise.sigma )

