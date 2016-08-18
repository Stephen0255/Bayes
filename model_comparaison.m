clear all

% Simulate the data

model_experiment.type = 'first_order_homogenous'; % first_order_homogenous, first_order, second order, Navier-Stokes
model_experiment.coefficients = [1, 2]; % 1, 2, 3, ... coefficients
coefficients=[1, 0.1, 0, 1, 2]; %  [y_initial, delta_t, t(0), t(end)]
noise.isnoise = true;
noise.sigma_noise = 0.05;


q_experiment=makedist('Normal','mu',0,'sigma',noise.sigma_noise); % distribution normale centr�e en 0

measurements = simulate_data(model_experiment, coefficients, noise, q_experiment); % measurements.y, measurements.t, measurements.number

plot(measurements.t,measurements.y)
%% Computation of the model


model = 'first_order_homogenous';
lower_a = 0; % borne inf�rieure pour la prior uniforme sur a
upper_a = 2; % borne sup�rieure
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

N = 1000; % nombres de d'�chantillons
M = 10000; % nombre d'it�rations pour la cha�ne de Markov
sigma_algorithm=0.1;
q_algorithm=makedist('Normal','mu',0,'sigma',sigma_algorithm); % distribution normale centr�e en 0

parameters.scheme = randi([-1 1],1,N); % upwind 1, downwind -1, center 0
parameters.several_scheme = true;
parameters.coefficients = zeros(number_coefficients,N); % premiere ligne: a, deuxieme: b, ...
parameters.poids = zeros(1,N);

for k=1:N % initialisation des valeurs
    
    if strcmp(model,'first_order_homogenous') % y' +ay = 0, un seul param�tre a
        parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
        y_calc=calcul_position(measurements,model,parameters.scheme(k), [parameters.coefficients(1,k)]); % calcul la position de tous les points
        
    elseif strcmp(model,'first_order') % y' + ay + b = 0, deux param�tres a et b
        parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
        parameters.coefficients(2,k) = random(p_b); % initialisation de la valeur de b
        y_calc=calcul_position(measurements,model,parameters.scheme(k),[parameters.coefficients(1,k) parameters.coefficients(2,k)]); % calcul la position de tous les points
        
    end
    somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carr�s des erreurs
    parameters.poids(k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant  

end

% initialisation
center = find(~parameters.scheme); % indices correspondant au sch�ma centr�
upwind = find(~(1-parameters.scheme)); % idem pour l'amont
downwind = find(~(1+parameters.scheme)); % idem pour l'aval

total_acceptance = 0;
total_acceptance_scheme = 0;

for i=1:M % � chaque it�ration, on fait :

    N=length(parameters.coefficients(1,:));
    compteur = 0;
    compteur_scheme = 0;

    
    for k=1:N % pour chaque �chantillon, on fait :
        scheme_post = randi([-1 1],1,1); % on tire un nouveau sch�ma
        if ~(scheme_post == parameters.scheme(k)) % si le sch�ma est diff�rent du pr�c�dent, on choisit un indice dont on prend les param�tres
            if scheme_post == 1 % Upwind
                n = length(upwind);
                index_scheme = randi([1 n],1,1);
                index = upwind(index_scheme);
            elseif scheme_post == 0 % Center
                n = length(center);
                index_scheme = randi([1 n],1,1);
                index = center(index_scheme);
            elseif scheme_post == -1 % Downwind
                n = length(downwind);
                index_scheme = randi([1 n],1,1);
                index = downwind(index_scheme);
            end
            
            parameters_post.coefficients = parameters.coefficients(:,index); % on donne les coefficients associ� � cet �chantillon
            parameters_post.poids = parameters.poids(index); % ainsi que le poids
        
            probability_acceptance=min(1,parameters_post.poids/parameters.poids(k));
        
            p=rand;
            if p<probability_acceptance
                parameters.coefficients(:,k) = parameters_post.coefficients;
                parameters.poids(k) = parameters_post.poids;
                compteur_scheme = compteur_scheme+1;
            end
            
        end
        
        [ parameters, compteur ] = markov_iteration( measurements, model, parameters, M, sigma_algorithm, q_algorithm, k, compteur );
    end
    
    acceptance=compteur/N;
    acceptance_scheme=compteur_scheme/N;
    total_acceptance = total_acceptance + acceptance;
    total_acceptance_scheme = total_acceptance_scheme + acceptance_scheme;
    calcul_effectue = i/M;
    avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
    display(avancement);
    
    center = find(~parameters.scheme); % indices correspondant au sch�ma centr�
    upwind = find(~(1-parameters.scheme)); % idem pour l'amont
    downwind = find(~(1+parameters.scheme)); % idem pour l'aval
end

acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/M)},{'%'});
display(acceptance_ratio);
acceptance_ratio_scheme = strcat({'Acceptance ratio scheme: '},{num2str(100*total_acceptance_scheme/M)},{'%'});
display(acceptance_ratio_scheme);

coefficients_upwind = parameters.coefficients(:,upwind);
coefficients_downwind = parameters.coefficients(:,downwind);
coefficients_center = parameters.coefficients(:,center);

figure()
histogram(parameters.scheme)
title('Distribution of the schemes');
xlabel('Downwind, Center, Upwind');


if strcmp(model,'first_order_homogenous')
    a_upwind=coefficients_upwind(1,:);
    a_downwind=coefficients_downwind(1,:);
    a_center=coefficients_center(1,:);
    
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
    
    a_upwind=coefficients_upwind(1,:);
    a_downwind=coefficients_downwind(1,:);
    a_center=coefficients_center(1,:);
    b_upwind=coefficients_upwind(2,:);
    b_downwind=coefficients_downwind(2,:);
    b_center=coefficients_center(2,:);    
    
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


    
    
    