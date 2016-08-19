function [ parameters ] = initialization_MCMC( measurements, model, scheme, prior, number_samples, sigma_algorithm, several_scheme)
% Initialise les param�tres de l'algorithme de M�tropolis-Hasting
%   Outpout : la valeur des coefficients pour le d�but de la cha�ne de
%   Markov et le poids (tout est regroup� sous parameters)
%   Inpout : les donn�es exp�rimentales, le mod�le ainsi que le sch�ma � utiliser. La prior sur les
%   coefficients du mod�le ainsi que le sigma de l'algorithme


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

if several_scheme % adapte la fonction pour traiter le cas o� il y a plusieurs sch�mas

    for k=1:number_samples
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
else
    for k=1:number_samples
        if strcmp(model,'first_order_homogenous') % y' +ay = 0, un seul param�tre a
            parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
            y_calc=calcul_position(measurements,model,parameters.scheme, [parameters.coefficients(1,k)]); % calcul la position de tous les points
        
        elseif strcmp(model,'first_order') % y' + ay + b = 0, deux param�tres a et b
            parameters.coefficients(1,k) = random(p_a); % initialisation de la valeur de a
            parameters.coefficients(2,k) = random(p_b); % initialisation de la valeur de b
            y_calc=calcul_position(measurements,model,parameters.scheme,[parameters.coefficients(1,k) parameters.coefficients(2,k)]); % calcul la position de tous les points
        
        end
        somme_carree = sum((measurements.y-y_calc).^2); % calcul la somme des carr�s des erreurs
        parameters.poids(k)=1/(sigma_algorithm^(measurements.number))*exp(-1/(2*sigma_algorithm^2)*somme_carree); % donne la valeur de PI pour la a correspondant  
    end    
    
end

end

