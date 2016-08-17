function [ parameters, acceptance ] = markov_iteration( measurements, model, scheme, parameters, M, sigma_algorithm, q )
% This function proceed to one iteration of the Metropolis-Hasting
% algorithm
%   The input arguments are the data (measurements), the model
%   (first_order, second_order, Navier-Stokes), the scheme (upwind, center,
%   ...) and the set of parameters.
%   The outpout will be a new set of parameters

N=length(parameters(1,:));

compteur = 0;
    
for k=1:N % pour chaque échantillon, on fait :
    
    if strcmp(model,'first_order_homogenous')
        
        a_post=parameters(1,k)+random(q); % choisi un nombre aléatoire centré en a(k)
        
        y_calc_post=calcul_position(measurements,model,scheme,[a_post]); % calcul la position de tout les points pour la nouvelle valeur
        somme_carree_post = sum((measurements.y-y_calc_post).^2);
        
        PI_post=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree_post);
        
        parameters_post=[a_post; PI_post];
      
    elseif strcmp(model,'first_order')
        
        a_post=parameters(1,k)+random(q); % choisi un nombre aléatoire centré en a(k)
        b_post=parameters(2,k)+random(q); % choisi un nombre aléatoire centré en b(k)
        
        y_calc_post=calcul_position(measurements,model,scheme,[a_post b_post]); % calcul la position de tout les points pour la nouvelle valeur
        somme_carree_post = sum((measurements.y-y_calc_post).^2);
        
        PI_post=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree_post);
        
        parameters_post=[a_post; b_post; PI_post];
        
    end
        
    probability_acceptance=min(1,PI_post/parameters(end,k));
        
    p=rand;
    if p<probability_acceptance
        parameters(:,k) = parameters_post;
        compteur = compteur+1;
    end
end

acceptance=compteur/N;

end

