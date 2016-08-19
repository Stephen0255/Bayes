function [ parameters ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme)
% Initialise les paramètres de l'algorithme de Métropolis-Hasting
%   Outpout : la valeur des coefficients pour le début de la chaîne de
%   Markov et le poids (tout est regroupé sous parameters)
%   Inpout : les données expérimentales, le modèle ainsi que le schéma à utiliser. La prior sur les
%   coefficients du modèle ainsi que le sigma de l'algorithme


if strcmp(model,'first_order_homogenous')
    number_coefficients = 1;
    p_a = makedist('Uniform','lower',prior.lower_a,'upper',prior.upper_a); % distribution uniforme (prior) pour a
elseif strcmp(model,'first_order')
    number_coefficients = 2;
    p_a = makedist('Uniform','lower',prior.lower_a,'upper',prior.upper_a); % distribution uniforme (prior) pour a
    p_b = makedist('Uniform','lower',prior.lower_b,'upper',prior.upper_b); % distribution uniforme (prior) pour b
end

parameters.several_scheme = several_scheme;
parameters.coefficients = zeros(number_coefficients,number_samples); % premiere ligne: a, deuxieme: b, ...
parameters.poids = zeros(1,number_samples);
parameters.scheme = scheme;

if several_scheme % adapte la fonction pour traiter le cas où il y a plusieurs schémas

    for k=1:number_samples
        if strcmp(model,'first_order_homogenous') % y' +ay = 0, un seul paramètre a
            parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
            y_calc=calcul_position(measurements,model,parameters.scheme(k), [parameters.coefficients(1,k)]); % calcul la position de tous les points
        
        elseif strcmp(model,'first_order') % y' + ay + b = 0, deux paramètres a et b
            parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
            parameters.coefficients(2,k) = random(p_b); % initialisation de la valeur de b
            y_calc=calcul_position(measurements,model,parameters.scheme(k),[parameters.coefficients(1,k) parameters.coefficients(2,k)]); % calcul la position de tous les points
        
        end
        somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carrés des erreurs
        parameters.poids(k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant  
    end
else
    for k=1:number_samples
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
    
end

end

