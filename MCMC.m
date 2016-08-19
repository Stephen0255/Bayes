function [ parameters ] = MCMC( measurements, model, parameters, sigma_algorithm, scheme_comparaison, number_iteration )
% Proceed to the Metropolis-Hasting algorithm
%  Outpout : parameters.scheme, parameters.coefficients (si besoin,
%  parameters.coefficients_upwind, _downwind, _center)

q_algorithm=makedist('Normal','mu',0,'sigma',sigma_algorithm); % distribution normale centrée en 0

if ~scheme_comparaison
    
    total_acceptance = 0;

    for i=1:number_iteration % à chaque itération, on fait :
        N=length(parameters.coefficients(1,:));
        compteur = 0;
        
        for k=1:N % pour chaque échantillon, on fait :
            [ parameters, compteur ] = markov_iteration( measurements, model, parameters, number_iteration, sigma_algorithm, q_algorithm, k, compteur );
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
    
    % initialisation
    center = find(~parameters.scheme); % indices correspondant au schéma centré
    upwind = find(~(1-parameters.scheme)); % idem pour l'amont
    downwind = find(~(1+parameters.scheme)); % idem pour l'aval

    total_acceptance = 0;
    total_acceptance_scheme = 0;

    for i=1:number_iteration % à chaque itération, on fait :
        N=length(parameters.coefficients(1,:));
        compteur = 0;
        compteur_scheme = 0;
    
        for k=1:N % pour chaque échantillon, on fait :
            scheme_post = randi([-1 1],1,1); % on tire un nouveau schéma
            if ~(scheme_post == parameters.scheme(k)) % si le schéma est différent du précédent, on choisit un indice dont on prend les paramètres
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
            
                parameters_post.coefficients = parameters.coefficients(:,index); % on donne les coefficients associé à cet échantillon
                parameters_post.poids = parameters.poids(index); % ainsi que le poids
                probability_acceptance=min(1,parameters_post.poids/parameters.poids(k));
        
                p=rand;
                if p<probability_acceptance
                    parameters.coefficients(:,k) = parameters_post.coefficients;
                    parameters.poids(k) = parameters_post.poids;
                    compteur_scheme = compteur_scheme+1;
                end            
            end % fin du if pour le cas où l'on change de schéma
        
            % réalisation de l'itération
            [ parameters, compteur ] = markov_iteration( measurements, model, parameters, number_iteration, sigma_algorithm, q_algorithm, k, compteur );
        end
    
        acceptance=compteur/N;
        acceptance_scheme=compteur_scheme/N;
        total_acceptance = total_acceptance + acceptance;
        total_acceptance_scheme = total_acceptance_scheme + acceptance_scheme;
        calcul_effectue = i/number_iteration;
        avancement = strcat({'Avancement: '}, {num2str(100*calcul_effectue)}, {'%'});
        display(avancement);
    
        center = find(~parameters.scheme); % indices correspondant au schéma centré
        upwind = find(~(1-parameters.scheme)); % idem pour l'amont
        downwind = find(~(1+parameters.scheme)); % idem pour l'aval
    end

    acceptance_ratio = strcat({'Acceptance ratio: '},{num2str(100*total_acceptance/number_iteration)},{'%'});
    display(acceptance_ratio);
    acceptance_ratio_scheme = strcat({'Acceptance ratio scheme: '},{num2str(100*total_acceptance_scheme/number_iteration)},{'%'});
    display(acceptance_ratio_scheme);
    
    parameters.coefficients_upwind = parameters.coefficients(:,upwind);
    parameters.coefficients_downwind = parameters.coefficients(:,downwind);
    parameters.coefficients_center = parameters.coefficients(:,center);

end

end

