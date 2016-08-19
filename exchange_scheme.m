function [ parameters, compteur_scheme ] = exchange_scheme( parameters, sample, compteur_scheme, localization)
%   Change the scheme and choose new parameters according to the distribution
%   of the chosen scheme
%   Inpout : the parameters. Outpout : the new parameters

            scheme_post = randi([-1 1],1,1); % on tire un nouveau sch�ma
            if ~(scheme_post == parameters.scheme(sample)) % si le sch�ma est diff�rent du pr�c�dent, on choisit un indice dont on prend les param�tres
                if scheme_post == 1 % Upwind
                    n = length(localization.upwind);
                    index_scheme = randi([1 n],1,1);
                    index = localization.upwind(index_scheme);
                elseif scheme_post == 0 % Center
                    n = length(localization.center);
                    index_scheme = randi([1 n],1,1);
                    index = localization.center(index_scheme);
                elseif scheme_post == -1 % Downwind
                    n = length(localization.downwind);
                    index_scheme = randi([1 n],1,1);
                    index = localization.downwind(index_scheme);
                end
            
                parameters_post.coefficients = parameters.coefficients(:,index); % on donne les coefficients associ� � cet �chantillon
                parameters_post.poids = parameters.poids(index); % ainsi que le poids
                probability_acceptance=min(1,parameters_post.poids/parameters.poids(sample));
        
                p=rand;
                if p<probability_acceptance
                    parameters.scheme(sample) = scheme_post;
                    parameters.coefficients(:,sample) = parameters_post.coefficients;
                    parameters.poids(sample) = parameters_post.poids;
                    compteur_scheme = compteur_scheme+1;
                end            
            end % fin du if pour le cas o� l'on change de sch�ma


end

