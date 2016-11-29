%% Initialisation
clear all
close all


%% Param�tres des simulations
Plot_agents = 0; % Permet d'afficher la simulation (0 = pas afficher, 1 = afficher)

n_simulations = 50; % nombre de simulations
n_agents = 400; % nombre d'agents
n_time_step = 501; % nombre de steps de simulation
delta_t = 1/5; % pas de temps � chaque step de simulation

%% Param�tres des agents

N_informe_1 = 0;
poids_information_1 = 0;

distance_repulsion = 1;
distance_alignement = 6;
distance_attraction = 6;
angle_de_vue_total = 270;

vitesse = ones(1, n_agents) * 1; % vitesse constance
ecart_type_bruit = 5;% bruit sur le mouvement
angle_de_rotation_maximum = 2;


%% Matrice de r�sultats

MatriceEcartType = zeros(n_time_step, n_simulations);
MatriceMoyenne = zeros(n_time_step, n_simulations);
Alignement = zeros(n_time_step, n_simulations);
Precision = zeros(n_time_step, n_simulations);
%% Param�tres de l'espace
dimension_tore = 60
mat_ajustx = [zeros(1, n_agents) -ones(1,n_agents) zeros(1, n_agents) ones(1,n_agents) -ones(1,n_agents) ones(1,n_agents) -ones(1,n_agents) zeros(1, n_agents) ones(1,n_agents)]*dimension_tore;
mat_ajusty = [zeros(1, n_agents) ones(1,n_agents) ones(1,n_agents) ones(1,n_agents) zeros(1, n_agents) zeros(1, n_agents) -ones(1,n_agents) -ones(1,n_agents) -ones(1,n_agents)]*dimension_tore;

tic

angle_de_vue_demi = (angle_de_vue_total/2) * pi / 180;


direction_d_information = zeros(1, n_agents);

poids_information(1:N_informe_1) = poids_information_1;
poids_information(N_informe_1+1 : n_agents) = 0;

for s = 1 : 1 : n_simulations

    positionX_start = rand(1, n_agents) * dimension_tore; % position initiale de la tete des poissons (x)
    positionY_start = rand(1, n_agents) * dimension_tore; % position initiale de la tete des poissons (y)
    direction_start = rand(1, n_agents) * 2 * pi; % direction initiale poissons

    delta_new_direction = zeros(1, n_agents);

    positionX = repmat(positionX_start, 1, 9) + mat_ajustx;
    positionY = repmat(positionY_start, 1, 9) + mat_ajusty;
    direction = repmat(direction_start, 1, 9);


    for t = 1 : 1 : n_time_step

        for i = 1 : 1 : n_agents

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FISH

			%%% Translation pour centrer sur l agent i

            positionX_temporaire = positionX - positionX(i);
            positionY_temporaire = positionY - positionY(i);


			%%% Rotation sur centrer sur l agent i

            positionX_temporaire_rotation = positionX_temporaire * cos(direction(i)) + positionY_temporaire * sin(direction(i));
            positionY_temporaire_rotation = -positionX_temporaire * sin(direction(i)) + positionY_temporaire * cos(direction(i));


			%%% Direction des autres agent
            directionautre = mod(direction - direction(i), 2*pi);

            directionautre (directionautre>pi) = directionautre (directionautre>pi) - 2*pi;


			%%% Translation pour centrer sur l agent i

            [Angle_Agents, Distance_Agents] = cart2pol(positionX_temporaire_rotation, positionY_temporaire_rotation);

            Agents_Percus = ((mod(Angle_Agents,2*pi) < angle_de_vue_demi) & Distance_Agents < distance_attraction & Distance_Agents > 0) | ((mod(Angle_Agents,2*pi) > (2*pi - angle_de_vue_demi)) & Distance_Agents < distance_attraction & Distance_Agents > 0);
            Agents_Selected = find(Agents_Percus==1);


%%%%

            if numel(Agents_Selected) > 0
                if sum( Distance_Agents(Agents_Selected) < distance_repulsion) > 0
                    delta_new_direction(i) = circ_mean( ( atan2( positionY_temporaire_rotation(Agents_Selected(Distance_Agents(Agents_Selected) < distance_repulsion)), positionX_temporaire_rotation(Agents_Selected(Distance_Agents(Agents_Selected) < distance_repulsion) ) ) + pi)');
                else
                    delta_new_direction(i) = circ_mean([directionautre(Agents_Selected) atan2(positionY_temporaire_rotation(Agents_Selected) , positionX_temporaire_rotation(Agents_Selected))]');
                end
            else
                delta_new_direction(i) = 0;
            end

			%%%%%%% Information

		    delta_new_direction(i) = atan2(sin(delta_new_direction(i)) + sin(direction_d_information(i)-direction(i))*poids_information(i), cos(delta_new_direction(i)) + cos(direction_d_information(i)-direction(i))*poids_information(i));
%%%%

            %%%%%%%%%%%%% AFFICHAGE

            if Plot_agents==1 && i==1
                figure(1)
                quiver(positionX(1:n_agents), positionY(1:n_agents), vitesse.*cos(direction(1:n_agents)), vitesse .*sin(direction(1:n_agents)), 0)
                hold on
                plot(positionX(1:n_agents), positionY(1:n_agents), '.', 'MarkerSize',10)
                axis([0 dimension_tore 0 dimension_tore])
                axis square
                %set(gca, 'XTickLabel', '', 'YTickLabel', '')
                hold off
                pause(0.1)
            end
        end

        delta_new_direction(delta_new_direction>pi) = delta_new_direction(delta_new_direction>pi) - 2*pi;

        delta_new_direction(delta_new_direction > angle_de_rotation_maximum * delta_t) = angle_de_rotation_maximum * delta_t;
        delta_new_direction(delta_new_direction < -angle_de_rotation_maximum * delta_t) = -angle_de_rotation_maximum * delta_t;

        arand = randn(1,n_agents) * ecart_type_bruit;
        direction_new = direction(1:n_agents) + delta_new_direction + arand;

        %arand = rand(1,n_agents)*ecart_type_bruit - ecart_type_bruit/2;
        positionX_new = mod(positionX(1:n_agents) + delta_t * vitesse .* cos(direction_new(1:n_agents)) , dimension_tore);
        positionY_new = mod(positionY(1:n_agents) + delta_t * vitesse .* sin(direction_new(1:n_agents)) , dimension_tore);

        positionX = repmat(positionX_new, 1, 9) + mat_ajustx;
        positionY = repmat(positionY_new, 1, 9) + mat_ajusty;
        direction = repmat(direction_new, 1, 9);

        MatriceEcartType(t, s) = circ_std(direction(1:n_agents)');
        MatriceMoyenne(t, s) = circ_mean(direction(1:n_agents)');
        Alignement(t,s)= (1/n_agents)*sqrt(sum(cos(direction_new(1:n_agents)))^2 + sum(sin(direction_new(1:n_agents)))^2 );
        Precision(t, s)= 0.5 * sqrt((cos(0)+cos(circ_mean(direction(1:n_agents)')))^2 + (sin(0) + sin(circ_mean(direction(1:n_agents)')))^2);
    end

    disp(['Simulation ' num2str(s) ' sur ' num2str(n_simulations)]);



    toc
end

vecteur_temps = [delta_t: delta_t : delta_t*n_time_step];

fprintf(['Moyenne de l orientation  : ', num2str(mean( MatriceMoyenne(n_time_step, :) )),'\n'])
fprintf(['EcartType de l orientation : ', num2str(mean( MatriceEcartType(n_time_step, :) )),'\n'])
fprintf(['Alignement: ', num2str(mean( Alignement(n_time_step,:) )),'\n'])
fprintf(['Precision: ', num2str(mean( Precision(n_time_step,:) )),'\n'])


figure(2)
subplot(2,1,1)
plot(vecteur_temps, MatriceEcartType)
hold on
plot(vecteur_temps, mean(MatriceEcartType,2),'-r', 'linewidth',2)
xlabel( 'temps(s)' )%Titre de l'axe X
ylabel( 'EcartType de l orientation' )%Titre de l'axe Y
title( 'variation de l EcartType d orientation en fonction de temps')%Titre de figure
xlim([0 n_time_step * delta_t])

subplot(2,1,2)
plot(vecteur_temps, MatriceMoyenne)
hold on
plot(vecteur_temps, mean(MatriceMoyenne, 2),'-r', 'linewidth',3)
xlabel( 'temps(s)' )
ylabel( 'Moyenne de l orientation' )
title( 'variation de la Moyenne de l orientation en foncion de temps')
xlim([0 n_time_step * delta_t])

figure(3)
subplot(2,1,1)
plot(vecteur_temps, Alignement)
hold on
plot(vecteur_temps, mean(Alignement,2),'-r', 'linewidth',2)
xlabel( 'temps(s)' )%Titre de l'axe X
ylabel( 'Alignement' )%Titre de l'axe Y
title( 'variation d Alignement en fonction de temps')%Titre de figure
xlim([0 n_time_step * delta_t])

subplot(2,1,2)
plot(vecteur_temps, Precision)
hold on
plot(vecteur_temps, mean(Precision, 2),'-r', 'linewidth',3)
xlabel( 'temps(s)' )
ylabel( 'Precision' )
title( 'variation de la precision en foncion de temps')
xlim([0 n_time_step * delta_t])
