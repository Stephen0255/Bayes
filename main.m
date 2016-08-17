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

model_comparaison = false;
model = 'first_order_homogenous';
lower_a = 0; % borne inférieure pour la prior uniforme sur a
upper_a = 2; % borne supérieure
lower_b = 0; % idem pour b
upper_b = 4;


if strcmp(model,'first_order_homogenous')
    number_coefficients = 1;
    p_a = makedist('Uniform','lower',lower_a,'upper',upper_a); % distribution uniforme (prior) pour a
elseif strcmp(model,'first_order')
    number_coefficients = 2;
    p_a = makedist('Uniform','lower',lower_a,'upper',upper_a); % distribution uniforme (prior) pour a
    p_b = makedist('Uniform','lower',lower_b,'upper',upper_b); % distribution uniforme (prior) pour b
end

N = 250; % nombres de d'échantillons
M = 250; % nombre d'itérations pour la chaîne de Markov
sigma_algorithm=0.1;
q_algorithm=makedist('Normal','mu',0,'sigma',sigma_algorithm); % distribution normale centrée en 0

parameters.scheme = 1; % upwind 1, downwind -1, center 0
parameters.coefficients = zeros(number_coefficients,N); % premiere ligne: a, deuxieme: b, ...
parameters.poids = zeros(1,N);

for k=1:N
    if strcmp(model,'first_order_homogenous') % y' +ay = 0, un seul paramètre a
        parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
        y_calc=calcul_position(measurements,model,parameters.scheme, [parameters.coefficients(1,k)]); % calcul la position de tous les points
        
    elseif strcmp(model,'first_order') % y' + ay + b = 0, deux paramètres a et b
        parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
        parameters.coefficients(2,k) = random(p_b); % initialisation de la valeur de b
        y_calc=calcul_position(measurements,model,parameters.scheme,[parameters.coefficients(1,k) parameters.coefficients(2,k)]); % calcul la position de tous les points
        
    end
    somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carrés des erreurs
    parameters.poids(k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant  

end

total_acceptance = 0;

for i=1:M % à chaque itération, on fait :
    N=length(parameters.coefficients(1,:));

    compteur = 0;
    
    for k=1:N % pour chaque échantillon, on fait :
    [ parameters, compteur ] = markov_iteration( measurements, model, parameters, M, sigma_algorithm, q_algorithm, k, compteur );
    end
    
    acceptance=compteur/N;
    total_acceptance = total_acceptance + acceptance;
    calcul_effectue = i/M;
    avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
    display(avancement);
end

acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/M)},{'%'});
display(acceptance_ratio);




if strcmp(model,'first_order_homogenous')
    figure()
    histogram(parameters.coefficients(1,:),'BinWidth', noise.sigma_noise)
    title(parameters.scheme);
    xlabel('a');
elseif strcmp(model,'first_order')
    figure()
    histogram(parameters.coefficients(1,:),'BinWidth', noise.sigma_noise)
    title(parameters.scheme);
    xlabel('a');    
    figure()
    histogram(parameters.coefficients(2,:),'BinWidth', noise.sigma_noise)
    title(parameters.scheme);
    xlabel('b');
end

