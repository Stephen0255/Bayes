function [ y ] = calcul_position(measurements, model, scheme, parameters )
% Calcul la position de tous les y à partir de coefficient de l'équation
% différentielles et de la position initiale.
%   Renvoie y un tableau de longueur le nombre de mesures N.
%   Il faut donner les résultas expérimentaux
%   Il faut préciser le modèle: 1 ordre, 2 ordre, Navier-Stokes, Équation
%   de la chaleur, puis le schéma: upwind, downwind, centré et enfin les
%   paramètres de ce schéma (sous forme de tableau) [a b]

N = measurements.number;
delta_discretization = measurements.delta;
y = zeros(1,N);
y(1)=measurements.y_initial;


if strcmp(model,'first_order_homogenous') % y' +ay = 0
    a = parameters(1);
    if strcmp(scheme,'upwind')
        for k=2:N
            y(k)=(1-a*delta_discretization)*y(k-1);
        end
    elseif strcmp(scheme,'downwind')
        for k=2:N
            y(k)=y(k-1)/(1+a*delta_discretization);
        end
    elseif strcmp(scheme,'center')
        y(2)=y(1)*((1-a*delta_discretization)+1/(1+a*delta_discretization))/2; % besoin d'interpoler le premier terme
        for k=3:N
            y(k)=y(k-2)-2*parameters(1,k)*delta_discretization*y(k-1);
        end
    elseif strcmp(scheme,'center_1/2')
        for k=2:N
            y(k)=y(k-1)*((1-a*delta_discretization)+1/(1+a*delta_discretization))/2; % 1/2 upwind + 1/2 downwind
        end
    end       

elseif strcmp(model,'first_order') % y' +ay + b = 0
    a = parameters(1);
    b = parameters(2);
    if strcmp(scheme,'upwind')
        for k=2:N
            y(k)=(1-a*delta_discretization)*y(k-1)-b*delta_discretization;
        end
    elseif strcmp(scheme,'downwind')
        for k=2:N
            y(k)=(y(k-1)-b*delta_discretization)/(1+a*delta_discretization);
        end
    elseif strcmp(scheme,'center')
        y(2)=((1-a*delta_discretization)*y(1)-b*delta_discretization+(y(1)-b*delta_discretization)/(1+a*delta_discretization))/2; % besoin d'interpoler le premier terme
        for k=3:N
            y(k)=y(k-2)-2*a*delta_discretization*y(k-1)-2*b*delta_discretization;
        end
    elseif strcmp(scheme,'center_1/2')
        for k=2:N
            y(k)=((1-a*delta_discretization)*y(k-1)-b*delta_discretization+(y(k-1)-b*delta_discretization)/(1+a*delta_discretization))/2; % 1/2 upwind + 1/2 downwind
        end
    end        
end

end

