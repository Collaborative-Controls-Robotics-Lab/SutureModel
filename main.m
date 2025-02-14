function main()
clc
clear
warning('off', 'all');

% Authors: Kimia Forghani, Yancy Diaz-Mercado
%
% Date Created: 02/14/2025
%
% Copyright (c) 2025
% Collaborative Controls and Robotics Laboratory
% University of Maryland, College Park
%
% All rights reserved.
%init params
    timeStep = 0.048;           %Change this based on hardware
    n = 40;                     %Number of nodes in the chain
    ell_thread = n*5;           %Thread length
    deltaL=ell_thread/(n-1);
    needle_step=[0;0];
    separationEnforcementGain = 1e+1;
    SafetyDistance = 0.5;

%%
% Initial position of the needle
    needle_pos = [0; -40];          %[x; y]
    needle_pos_temp = needle_pos;

    initSep = 0.99*deltaL;          %Initial distanc between each two nodes
    initPos = needle_pos';
    threadXY = [ones(n, 1)*initPos(1,1)  ,  (initPos(1,2)-initSep: -initSep: (n-1)*-initSep+initPos(1,2)-initSep)'];
    threadXY=threadXY'; %size(threadXY) = 2 X n

    % Create a figure for capturing key presses
    hFig = figure('KeyPressFcn', @keyPress);

    hNeedle =  plot(needle_pos(1,1), needle_pos(2,1), 'r-o');
    hold on

    hThread = plot(threadXY(1,:), threadXY(2,:) , '-o','LineWidth',3,'MarkerSize', 2, 'MarkerFaceColor', 0.75*[1 1 1],'Color', 0.5*[1 1 1]); % Suture initially rest on x axis 
    axis([-80 80 -120 40 -5 5]);
    axis off; 
    set(gcf, 'Color', 'w');
    
%% Visualize Abdominal Wall and Inguinal Ring
%Visualization parameters
    rho_p = 0.00071/2*1000; %[mm] outer radius of the needle
    percentageFromEdge2Ring = 0.85;
%Hernia visualization
    [~,abdominalWall,~] = visualizeInguinalRing('NeedleRadius',2*rho_p,'FractionEdgeToRing',percentageFromEdge2Ring);

%Obstacle Triangulation
    tri_vertices{1} = triangles_vertices_delaunay(abdominalWall{1}); %project point method
    tri_vertices{2} = triangles_vertices_delaunay(abdominalWall{2}); %project point method
    tri_vertices{3} = triangles_vertices_delaunay(abdominalWall{3}); %project point method

%%
%Choose your barrier function

%Without Stiffness Lyapunov Function: suitable for silk suture
    barrierCertificate = create_si_connectivity_barrier_certificate_with_obstacles('MaxSeparation',deltaL,'SafetyDistance', SafetyDistance,'BarrierGain',separationEnforcementGain,'N',n,'tri_verices',tri_vertices);

%With Stiffness Lyapunov Function: suitable for Polyamide suture
    %barrierCertificate = create_si_connectivity_barrier_certificate_with_obstacles_stiff('MaxSeparation',deltaL,'SafetyDistance', SafetyDistance,'BarrierGain',separationEnforcementGain,'N',n,'tri_verices',tri_vertices);

%%
%initial thread vlocity
    threadXYVel = zeros(2,n); %Initially stationary 
j=1;
tic; %start timer
while j<5000

    threadXYVel = 0.9* threadXYVel; %Damping effect Beta=0.9
    du = (needle_step)./timeStep;   %Needle velocity based on user input
    needle_step =[0;0];             %Change user input back to zero
   
%************************ Thread Update Start *************************
    [threadXYVel, duC, B(j,:), Btis(j,:)] = barrierCertificate(threadXYVel,threadXY,needle_pos_temp,du);

%CLF 
    if j>1
        % slows down needle if two consecutive obs
        Btis(j,:) = Btis(j,:)*1e-3;
    
        % Define the lengths of the sections
        n1 = n+1;
    
        % Extract the three parts of the vector
        Btis_part1 = Btis(j, 1:n1);                   % First part
        Btis_part2 = Btis(j, n1+1:2*n1);              % Second part
        Btis_part3 = Btis(j, 2*n1+1:3*n1);            % Third part
    
        Bthresh = 0.03;
        velthresh = 4;
        % Check if Btis(j, i) in any two of the three parts is < Bthresh
        condition1 = (Btis_part1 < Bthresh & Btis_part2 < Bthresh) | ...
                     (Btis_part1 < Bthresh & Btis_part3 < Bthresh) | ...
                     (Btis_part2 < Bthresh & Btis_part3 < Bthresh);
                 if any(condition1)  
                        % Find the index where the condition is true
                        index = find(condition1, 2); % Find the first index where condition1 is true
        
                        % Check if condition1 is true for any element
                        if ~isempty(index) && any(arrayfun(@(i) norm(threadXYVel(:, mod(i,n)+1)), index) > velthresh) 
                            duC = duC * 0.9^sum(condition1);
                            threadXYVel = threadXYVel* 0.9^sum(condition1);
                            hThread.Color = [1, 165/255 * (13 - sum(condition1)) / 12,0]; %orange
                        else 
                            hThread.Color = [0,0.2,1]; %blue
                        end
        
                 else
                         hThread.Color = [0,1,0]; %green
                 end
    end

    %Update thread visual
    [threadXY, threadXYVel, hThread ,needle_pos_temp] = ThreadUpdateBarrier(threadXY, threadXYVel, timeStep, hThread,needle_pos_temp, duC);
        
    %************************* Thread Update End **************************

    %Update needle 
    hNeedle.XData = needle_pos_temp(1,1);
    hNeedle.YData = needle_pos_temp(2,1);
    
    % Update the figure window 
    drawnow;  
%%
         elapsedTime = toc;  % Check the elapsed time
       
    if elapsedTime < timeStep
        pause(timeStep - elapsedTime);  % Adjust the pause to maintain consistent timing
    else
        disp(elapsedTime)
    end
  tic;  % Reset timer after each loop iteration

j=j+1;
end    

 figure()
    plot(1:j-1,B(:,:))
    title('h con1')
    xlabel('Iteration')
    ylabel('h(xi)')
%     legend('node 1','node 10','node 20','node 30','node 40','node 50')


%     figure()
%     plot(2:j-1,1e3*Btis(2:j-1,:))
%     xlabel('Iteration')
%     ylabel('h(xi)')
%     title('h obs1')

    figure()
    plot(2:j-1,1e3*Btis(2:j-1,1:n+1))
    xlabel('Iteration')
    ylabel('h(xi)')
    title('h obs1')

    figure()
    plot(2:j-1,1e3*Btis(2:j-1,n+2:2*n+2))
    xlabel('Iteration')
    ylabel('h(xi)')
    title('h obs2')

        figure()
    plot(2:j-1,1e3*Btis(2:j-1,2*n+3:3*n+3))
    xlabel('Iteration')
    ylabel('h(xi)')
    title('h obs3')

function keyPress(~, event)
        step_size = 0.475;  % How much the needle moves with each key press 
        switch event.Key
            case 'w'  % Move needle up
                needle_step = + step_size.*[0;1];
            case 's'  % Move needle down
                needle_step = - step_size.*[0;1];
            case 'a'  % Move needle left
                needle_step = - step_size.*[1;0];
            case 'd'  % Move needle right
                needle_step = + step_size.*[1;0];
                case 'q'  % Move needle right
                needle_step = + step_size.*[0;0];
        end

    end 
end
