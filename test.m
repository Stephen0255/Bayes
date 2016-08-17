clear all

% Simulate the data

model_experiment.type = 'first_order_homogenous'; % first_order_homogenous, first_order, second order, Navier-Stokes
model_experiment.coefficients = [1, 2]; % 1, 2, 3, ... coefficients
coefficients=[1, 0.1, 0, 1, 2]; %  [y_initial, delta_t, t(0), t(end)]
noise.isnoise = true;
noise.sigma_noise = 0.01;


q=makedist('Normal','mu',0,'sigma',noise.sigma_noise); % distribution normale centr�e en 0

measurements = simulate_data(model_experiment, coefficients, noise, q); % measurements.y, measurements.t, measurements.number

%


model = 'first_order_homogenous';
scheme='upwind'; % upwind, downwind, center
N = 1000; % nombres de d'�chantillons
M = 1000; % nombre d'it�rations pour la cha�ne de Markov
sigma_algorithm=noise.sigma_noise;

if strcmp(model,'first_order_homogenous') % y' +ay = 0, un seul param�tre a
    p_a = makedist('Uniform','lower',0,'upper',2); % distribution uniforme (prior) pour a
    parameters = zeros(2,N); % premi�re ligne: a, deuxi�me ligne, le poids de a: PI
    
    for k=1:N
    parameters(1,k) = random(p_a); % initialisation de la valeur de a
    
    y_calc=calcul_position(measurements,model,scheme, [parameters(1,k)]); % calcul la position de tous les points
    somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carr�s des erreurs
    
    parameters(2,k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant  
    end
    
elseif strcmp(model,'first_order') % y' + ay + b = 0, deux param�tres a et b    
    p_a = makedist('Uniform','lower',0,'upper',2); % distribution uniforme (prior) pour a
    p_b = makedist('Uniform','lower',0,'upper',4); % distribution uniforme (prior) pour b
    parameters = zeros(3,N); % premi�re ligne a, deuxi�me ligne b, troisi�me ligne PI
    
    for k=1:N
    parameters(1,k) = random(p_a); % initialisation de la valeur de a
    parameters(2,k) = random(p_b); % initialisation de la valeur de b
    
    y_calc=calcul_position(measurements,model,scheme,[parameters(1,k) parameters(2,k)]); % calcul la position de tous les points
    somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carr�s des erreurs
    
    parameters(3,k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant      
    end    
end

total_acceptance = 0;

for i=1:M % � chaque it�ration, on fait :
[ parameters, acceptance ] = markov_iteration( measurements, model, scheme, parameters, M, sigma_algorithm, q );

total_acceptance = total_acceptance + acceptance;
calcul_effectue = i/M;
avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
display(avancement);
end

acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/M)},{'%'});
display(acceptance_ratio);




if strcmp(model,'first_order_homogenous')
    figure()
    histogram(parameters(1,:))
    title(scheme);
    xlabel('a');
elseif strcmp(model,'first_order')
    figure()
    histogram(parameters(1,:))
    title(scheme);
    xlabel('a');    
    figure()
    histogram(parameters(2,:))
    title(scheme);
    xlabel('b');
end

