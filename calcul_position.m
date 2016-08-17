function [ y ] = calcul_position(measurements, model, schema, parameters )
% Calcul la position de tout les y � partir de coefficient de l'�quation
% diff�rentielles et de la position initiale.
%   Renvoie y un tableau de longueur le nombre de mesures N.
%   Il faut donner les r�sultas exp�rimentaux
%   Il faut pr�ciser le mod�le: 1 ordre, 2 ordre, Navier-Stokes, �quation
%   de la chaleur, puis le sch�ma: upwind, downwind, centr� et enfin les
%   param�tres de ce sch�ma

N = measurements.number;
delta_discretization = measurements.delta;
y = zeros(1,N);
y(1)=measurements.y_initial;

if strcmp(model,'first_order')
    a = parameters(1); % y' +ay + b = 0
    b = parameters(2);
    if strcmp(schema,'upwind')
        for k=2:N
            y(k)=(1-a*delta_discretization)*y(k-1)-b*delta_discretization;
        end
    elseif strcmp(schema,'downwind')
        for k=2:N
            y(k)=(y(k-1)-b*delta_discretization)/(1+a*delta_discretization);
        end
    elseif strcmp(schema,'center')
        y(2)=((1-a*delta_discretization)*y(1)-b*delta_discretization+(y(1)-b*delta_discretization)/(1+a*delta_discretization))/2; % besoin d'interpoler le premier terme
        for k=3:N
            y(k)=y(k-2)-2*a*delta_discretization*y(k-1)-2*b*delta_discretization;
        end
    elseif strcmp(schema,'center_1/2')
        for k=2:N
            y(k)=((1-a*delta_discretization)*y(k-1)-b*delta_discretization+(y(k-1)-b*delta_discretization)/(1+a*delta_discretization))/2; % 1/2 upwind + 1/2 downwind
        end
    end        
end

end

