function [ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, scheme_comparaison, number_iteration )
% Proceed to the Metropolis-Hasting algorithm
%  Outpout : parameters.scheme, parameters.coefficients (si besoin,
%  parameters.coefficients_upwind, _downwind, _center)

q_algorithm=makedist('Normal','mu',0,'sigma',sigma_algorithm); % distribution normale centrée en 0

if ~scheme_comparaison
    
    total_acceptance = 0;
    N=length(parameters.coefficients(1,:));

    for i=1:number_iteration % à chaque itération, on fait :
        compteur = 0;
        
        for k=1:N % pour chaque échantillon, on fait :
            [ parameters, compteur ] = markov_iteration( measurements, model, parameters, sigma_algorithm, q_algorithm, k, compteur );
        end 
        acceptance=compteur/N;
        total_acceptance = total_acceptance + acceptance;
        calcul_effectue = i/number_iteration;
        avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
        display(avancement);
    end

    acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/number_iteration)},{'%'});
    display(acceptance_ratio);
    
else
    
    % Localize each scheme
    localization.center = find(~parameters.scheme); % indices correspondant au schéma centré
    localization.upwind = find(~(1-parameters.scheme)); % idem pour l'amont
    localization.downwind = find(~(1+parameters.scheme)); % idem pour l'aval
    
    N=length(parameters.coefficients(1,:));

    total_acceptance = 0;
    total_acceptance_scheme = 0;

    for i=1:number_iteration % à chaque itération, on fait :
        compteur = 0;
        compteur_scheme = 0;
    
        for k=1:N % pour chaque échantillon, on fait :
            
            [ parameters, compteur_scheme ] = exchange_scheme( parameters, k, compteur_scheme, localization);
        
            % réalisation de l'itération
            [ parameters, compteur ] = markov_iteration( measurements, model, parameters, sigma_algorithm, q_algorithm, k, compteur );
        end
    
        acceptance=compteur/N;
        acceptance_scheme=compteur_scheme/N;
        total_acceptance = total_acceptance + acceptance;
        total_acceptance_scheme = total_acceptance_scheme + acceptance_scheme;
        calcul_effectue = i/number_iteration;
        avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
        display(avancement);
    
        localization.center = find(~parameters.scheme); % indices correspondant au schéma centré
        localization.upwind = find(~(1-parameters.scheme)); % idem pour l'amont
        localization.downwind = find(~(1+parameters.scheme)); % idem pour l'aval
    end

    acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/number_iteration)},{'%'});
    display(acceptance_ratio);
    acceptance_ratio_scheme = strcat({'Acceptance ratio scheme: '},{num2str(100*total_acceptance_scheme/number_iteration)},{'%'});
    display(acceptance_ratio_scheme);
    
    parameters.coefficients_upwind = parameters.coefficients(:,localization.upwind);
    parameters.coefficients_downwind = parameters.coefficients(:,localization.downwind);
    parameters.coefficients_center = parameters.coefficients(:,localization.center);

end

end

