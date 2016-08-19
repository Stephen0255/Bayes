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
choosen_scheme = 'Upwind'; % Upwind, Center or Downwind

prior.lower_a = 0; % borne inférieure pour la prior uniforme sur a
prior.upper_a = 2; % borne supérieure
prior.lower_b = 0; % idem pour b
prior.upper_b = 4;

if  ~scheme_comparaison
    if strcmp(choosen_scheme,'Upwind')
        scheme = 1; % upwind 1, downwind -1, center 0
    elseif strcmp(choosen_scheme,'Center')
        scheme = 0;
    elseif strcmp(choosen_scheme,'Downwind')
        scheme = -1;
    end
    several_scheme = false;
    
    [ parameters ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);
    [ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, several_scheme, number_iteration );
    create_histogram( model, parameters, several_scheme, noise.sigma )
    
else
    several_scheme = false;
    scheme = 1;
    [ parameters_upwind ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);
    [ parameters_upwind ] = MCMC( measurements, model, parameters_upwind, sigma_algorithm, several_scheme, number_iteration );
    
    scheme = 0;
    [ parameters_center ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);
    [ parameters_center ] = MCMC( measurements, model, parameters_center, sigma_algorithm, several_scheme, number_iteration );
    
    scheme = -1;
    [ parameters_downwind ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme);
    [ parameters_downwind ] = MCMC( measurements, model, parameters_downwind, sigma_algorithm, several_scheme, number_iteration );
    
    several_scheme = true;
    parameters.several_scheme = several_scheme;
    parameters.coefficients = [parameters_upwind.coefficients, parameters_center.coefficients, parameters_downwind.coefficients];
    parameters.poids = [parameters_upwind.poids, parameters_center.poids, parameters_downwind.poids];
    parameters.scheme = [ones(1,number_samples), zeros(1,number_samples), -1*ones(1,number_samples)];
    [ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, several_scheme, number_iteration );
end


create_histogram( model, parameters, several_scheme, noise.sigma )



