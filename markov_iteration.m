function [ parameters, compteur ] = markov_iteration( measurements, model, parameters, M, sigma_algorithm, q, iteration, compteur )
% This function proceed to one iteration of the Metropolis-Hasting
% algorithm
%   The input arguments are the data (measurements), the model
%   (first_order, second_order, Navier-Stokes), the scheme (upwind, center,
%   ...) and the set of parameters.
%   The outpout will be a new set of parameters

if model_comparaison
    scheme = parameters.scheme(k);
else
    scheme = parameters.scheme;
end

    
    if strcmp(model,'first_order_homogenous')
        
        a_post=parameters.coefficients(1,iteration)+random(q); % choisi un nombre aléatoire centré en a(k)
        parameters_post.coefficients = [a_post];
        
        y_calc_post=calcul_position(measurements,model,scheme,parameters_post.coefficients); % calcul la position de tout les points pour la nouvelle valeur
        somme_carree_post = sum((measurements.y-y_calc_post).^2);
        
        parameters_post.poids=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree_post);
      
    elseif strcmp(model,'first_order')
        
        a_post=parameters.coefficients(1,iteration)+random(q); % choisi un nombre aléatoire centré en a(k)
        b_post=parameters.coefficients(2,iteration)+random(q); % choisi un nombre aléatoire centré en b(k)
        
        parameters_post.coefficients = [a_post; b_post];
        
        y_calc_post=calcul_position(measurements,model,scheme,parameters_post.coefficients); % calcul la position de tout les points pour la nouvelle valeur
        somme_carree_post = sum((measurements.y-y_calc_post).^2);
        
        parameters_post.poids=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree_post);
        
    end
        
    probability_acceptance=min(1,parameters_post.poids/parameters.poids(iteration));
        
    p=rand;
    if p<probability_acceptance
        parameters.coefficients(:,iteration) = parameters_post.coefficients;
        parameters.poids(iteration) = parameters_post.poids;
        compteur = compteur+1;
    end

end

