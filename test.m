% Simule les donn�es pour y' + ay + b = 0

a=1;
b=2;
measurements.y_initial = 1;
measurements.delta=0.1;
t= 0:measurements.delta:1;
y=(measurements.y_initial+b/a)*exp(-a*t)-b/a;
plot(t,y)

% Rajouter une perturbation gaussienne

measurements.y=y;
measurements.t=t;
measurements.sigma=0.1;

clearvars y t

% Supposer une r�partition pour a

N = 500; % nombres de d'�chantillons
M = 10000; % nombre d'it�rations pour la cha�ne de Markov
measurements.number = length(measurements.y);
p_a = makedist('Uniform','lower',0,'upper',2); % distribution uniforme (prior) pour a
p_b = makedist('Uniform','lower',0,'upper',4); % distribution uniforme (prior) pour b
a = zeros(1,N); % valeurs de mon coefficient
b = zeros(1,N); % valeurs de mon coefficient
PI=zeros(1,N); % poids de l'erreur li� au coefficient
n = 2; % coefficient d'augmentation de la discr�tisation en temps (pour �tre plus pr�cis)
delta = measurements.delta/n;
sigma=measurements.sigma;

q=makedist('Normal','mu',0,'sigma',measurements.sigma);

% Initialisation de la valeur de chaque a
for k=1:N
    a(k) = random(p_a);
    b(k) = random(p_b);
end

model='first_order';
schema='upwind'; % upwind, downwind, center, center_1/2

tic
for i=1:M % � chaque it�ration, on fait :
    compteur = 0;
    
    for k=1:N % pour chaque �chantillon, on fait :
        a_post=a(k)+random(q); % choisi un nombre al�atoire centr� en a(k)
        b_post=b(k)+random(q); % choisi un nombre al�atoire centr� en b(k)
        
        y_calc=calcul_position(measurements,model,schema,[a(k), b(k)]); % calcul la position de tout les points pour l'ancienne valeur
        y_calc_post=calcul_position(measurements,model,schema,[a_post, b_post]); % calcul la position de tout les points pour la nouvelle valeur
        
        somme_carree = sum((measurements.y-y_calc).^2);
        somme_carree_post = sum((measurements.y-y_calc_post).^2);
        
        PI(k)=1/(sigma^(measurements.number))*exp(-1/(2*sigma^2)*somme_carree);
        PI_post=1/(sigma^(measurements.number))*exp(-1/(2*sigma^2)*somme_carree_post);
        
        acceptance_ratio=min(1,PI_post/PI(k));
        
        p=rand;
        if p<acceptance_ratio
            a(k)=a_post;
            b(k)=b_post;
            compteur = compteur+1;
        end
    end
    ratio_acceptance=compteur/N
    calcul_effectue=i/M
    
end
temps_markov = toc

figure()
histogram(a)
title(schema);
xlabel('a');

figure
histogram(b)
title(schema);
xlabel('b');





