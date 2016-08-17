function [ y ] = calcul_position(measurements, model, scheme, parameters )
% Calcul la position de tous les y � partir de coefficient de l'�quation
% diff�rentielles et de la position initiale.
%   Renvoie y un tableau de longueur le nombre de mesures N.
%   Il faut donner les r�sultas exp�rimentaux
%   Il faut pr�ciser le mod�le: 1 ordre, 2 ordre, Navier-Stokes, �quation
%   de la chaleur, puis le sch�ma: upwind, downwind, centr� et enfin les
%   param�tres de ce sch�ma (sous forme de tableau) [a b]

N=measurements.number;
delta=measurements.delta;

y = zeros(1,N);
y(1)=measurements.y_initial;


if strcmp(model,'first_order_homogenous') % y' +ay = 0
    a = parameters(1);
    if scheme == 1 % Upwind
        for k=2:N
            y(k)=(1-a*delta)*y(k-1);
        end
    elseif scheme == -1 % Downwind
        for k=2:N
            y(k)=y(k-1)/(1+a*delta);
        end
    elseif scheme == 0 % Center
        for k=2:N
            y(k)=y(k-1)*((1-a*delta)+1/(1+a*delta))/2; % 1/2 upwind + 1/2 downwind
        end
    end       

elseif strcmp(model,'first_order') % y' +ay + b = 0
    a = parameters(1);
    b = parameters(2);
    if scheme == 1 % Upwind
        for k=2:N
            y(k)=(1-a*delta)*y(k-1)-b*delta;
        end
    elseif scheme == -1 % Downwind
        for k=2:N
            y(k)=(y(k-1)-b*delta)/(1+a*delta);
        end
    elseif scheme == 0 % Center
        for k=2:N
            y(k)=((1-a*delta)*y(k-1)-b*delta+(y(k-1)-b*delta)/(1+a*delta))/2; % 1/2 upwind + 1/2 downwind
        end
    end        
end

end

