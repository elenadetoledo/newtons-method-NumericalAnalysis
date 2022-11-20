clear
clc
% NUMERICAL METHODS (PRACTICE 2)

% Miguel González Martínez (NIA: 100451423)
% Ricardo Macías Leiva (NIA: 100452215)
% Elena de Toledo Hernández (NIA: 100452170)

% PROGRAM 

% Write a program that implements the Newton-Raphson method in order to find
% zeros of the gradient of the Himmelblau’s function from a given starting 
% point x ∈ [−6, 6]×[−6, 6].

fprintf('\n\nPRACTICE 2: MULTIVARIABLE NEWTON_RAPHSON METHOD IN HIMMELBLAUs FUNCTION\n\n')

% The 3 points that were given to our team (Team 11)
P1 = [-6 6]; 
P2 = [0.01 -3.01];
P3 = [3.02 -3.03];

%GRID SEARCH -- Plotting the function
[x1, y1] = meshgrid(-6:0.01:6, -6:0.01:6);
Himm = himmelblaus(x1,y1); %Evaluate function at x1,y1
contour(x1,y1,Himm,100)
tolerance = 10^-6; % we use this value for tolerance because it has a good balance between effectiveness and efficiency. We use 6 decimals

% We need maxima, minima and saddle points of Himmelblau's to be able to plot the function. 
% Maxima and minimas are found in Wikipedia, Google 
% Saddle points found on the internet too

min1 = [3,2];
min2 = [-2.805118, 3.131312];
min3 = [-3.779310, -3.283186];
min4 = [3.584428, -1.848126];
max = [-0.270845, -0.923039];
sad1 =[-3.073030 -0.081353];
sad2 = [-0.127961 -1.953710];
sad3 = [0.086677 2.884250];
sad4 = [3.38515 0.073851];

%Create 2 loops so that every point contained in the rectangle
%[-6,6]x[-6,6] is taken into account
for i = -6:0.1:6
    for j = -6:0.1:6
        initialp = [i , j];
        xinitial = initialp(1); %%setting the coordinates of the initial point
        yinitial = initialp(2);
        error = 1;  %%use any value larger than tolerance so that the loop is entered
        endLoop = 0; %condition to access the while loop
        its = 0; %start with 0 its
        
        while (error > tolerance) && (endLoop == 0)
            its = its + 1; %increase number of its
            currX  = newtonraphson(xinitial, yinitial); %use newtonraphson
            error = norm(currX - initialp); %compute error
            initialp(1) = currX(1); %update values
            initialp(2) = currX(2);
            
            if its >= 10000 %use a maximal number of its so that the loop ends at some point in case of divergence
                exit = 1;
            end
        end
        %STUDY OF CONVERGENC. For this we study the
        %distance between between each of the points we got through Newton
        %method 
        %and the critical ones.
        
        dist1 = norm(currX-min1);
        dist2 = norm(currX-min2);
        dist3 = norm(currX-min3);
        dist4 = norm(currX-min4);
        dist5 = norm(currX-max);
        dist6 = norm(currX-sad1);
        dist7 = norm(currX-sad2);
        dist8 = norm(currX-sad3);
        dist9 = norm(currX-sad4);
        
        %We now must store all the distances computed above in a vector.
        vector = [dist1, dist2, dist3, dist4, dist5, dist6, dist7, dist8, dist9];
        
        % Now that all distances are stored, we compute the minimum one and plot it
        
        if min(vector) == dist1
                hold on
                plot(i,j, '--.k')
        elseif min(vector) == dist2
                hold on
                plot(i,j, '--.b')
        elseif min(vector) == dist3
                hold on
                plot(i,j, '--.m')
        elseif min(vector) == dist4
                hold on
                plot(i,j, '--.c')
        elseif min(vector) == dist5
                hold on
                plot(i,j, '--.y')
        elseif (min(vector) == dist6) || (min(vector) == dist7) || (min(vector) == dist8) || (min(vector) == dist9)
                hold on
                plot(i,j,'--.r')
        end
        vector = []; %Vector is emptied because it will be filled with other data in next iteration
    end
end

%We include the following figures in the plot so that it is understood
%better.
%On the minima, you will find an x
hold on
plot(min1(1), min1(2), '--xk')
hold on
plot(min2(1), min2(2), '--xk')
hold on
plot(min3(1), min3(2), '--xk')
hold on
plot(min4(1), min4(2), '--xk')
hold on
%On the maxima, you will find a circle
plot(max(1), max(2), '--ok')
hold on
%On the saddle points, you will see a cross 
plot(sad1(1), sad1(2), '--+k')
hold on
plot(sad2(1), sad2(2), '--+k')
hold on
plot(sad3(1), sad3(2), '--+k')
hold on
plot(sad4(1), sad4(2), '--+k')
hold off

%Now we should study in depth the three point P1, P2 and P3 that were given
%to us: 


fprintf('For the first point, P1: (-6, 6)\n');   
its = 0;     %Initialize the its to 0
fin= 0;    %condition to enter the while loop
error = 1;   %Choose any number greater than tolerance so that we may enter the loop
inpoint = P1; %set point X0 to P1
x0 = inpoint(1);  %Coordinate x of point
y0 = inpoint(2);  %Coordinate y of point

 %Two conditions must be satisfied for the while loop to end --> the error
 %is less than the tolerance or the exit condition is satisfied. This
 %will happen when the norm of the vector of distance between two
 %consecutive points is greater than the tolerance,
while (error > tolerance) && (fin == 0)
 J = jacobian(x0,y0); %for the condition number we need the Jacobian --> use the function we have created for said purpose   
 conditional = cond(J);  %Compute the conditional number of the jacobian matrix in each point          
its = its + 1; %%increase by 1 the number of its   
currX = newtonraphson(x0,y0); %Use the external function we have created to use the newton raphson and calculate following point
error = norm(currX-inpoint);  %compute the error
inpoint = currX;  %Update for Newton Raphson
x0 = inpoint(1);
y0 = inpoint(2);
%%use a maximal number of its so that the loop ends at some point in case of divergence
    if its >= 50
        fin = 1;
        disp('The point (-6,6) diverges') 
    end
end
fprintf('Condition number in last iteration was %.6f\n',conditional);
fprintf('Number of its needed: %d\n', its);
%If exit is equal to zero, at some point in the loop error <= tolerance -->
%convergence 
if fin == 0
   disp('The point (-6,6) converges')
end

fprintf('For the second point, P2: (0.01, 3.01)\n');  
its = 0;     %Initialize the its to 0
fin= 0;    %condition to enter the while loop
error = 1;   %Choose any number greater than tolerance so that we may enter the loop
inpoint = P2; %set point X0 to P1
x0 = inpoint(1);  %Coordinate x of point
y0 = inpoint(2);  %Coordinate y of point

 %Two conditions must be satisfied for the while loop to end --> the error
 %is less than the tolerance or the exit condition is satisfied. This
 %will happen when the norm of the vector of distance between two
 %consecutive points is greater than the tolerance,
while (error > tolerance) && (fin == 0)
 J = jacobian(x0,y0); %for the condition number we need the Jacobian --> use the function we have created for said purpose   
 conditional = cond(J);  %Compute the conditional number of the jacobian matrix in each point          
its = its + 1; %%increase by 1 the number of its   
currX = newtonraphson(x0,y0); %Use the external function we have created to use the newton raphson and calculate following point
error = norm(currX-inpoint);  %compute the error
inpoint = currX;  %Update 
x0 = inpoint(1);
y0 = inpoint(2);
%%use a maximal number of its so that the loop ends at some point in case of divergence
    if its >= 50
        fin = 1;
        disp('The point (0.01,3.01) diverges') 
    end
end
fprintf('Condition number in last iteration was %.6f\n',conditional);
fprintf('Number of its needed: %d\n', its);
%Method converges:
if fin == 0
   disp('The point (0.01,3.01) converges')
end
        
fprintf('For the third point, P3: (3.02 , -3.03)\n');  
its = 0;     %Initialize the its to 0
fin= 0;    %condition to enter the while loop
error = 1;   %Choose any number greater than tolerance so that we may enter the loop
inpoint = P3; %set point X0 to P1
x0 = inpoint(1);  %Coordinate x of point
y0 = inpoint(2);  %Coordinate y of point

 %Two conditions must be satisfied for the while loop to end --> the error
 %is less than the tolerance or the exit condition is satisfied. This
 %will happen when the norm of the vector of distance between two
 %consecutive points is greater than the tolerance,
while (error > tolerance) && (fin == 0)
 J = jacobian(x0,y0); %for the condition number we need the Jacobian --> use the function we have created for said purpose   
 conditional = cond(J);  %Compute the conditional number of the jacobian matrix in each point          
its = its + 1; %%increase by 1 the number of its   
currX = newtonraphson(x0,y0); %Use the external function we have created to use the newton raphson and calculate following point
error = norm(currX-inpoint);  %compute the error
inpoint = currX;  %Update 
x0 = inpoint(1);
y0 = inpoint(2);
%%use a maximal number of its so that the loop ends at some point in case of divergence
    if its >= 50
        fin = 1;
        disp('The point (3.02,-3.03) diverges') 
    end

end
fprintf('Condition number in last iteration was %.6f\n',conditional);
fprintf('Number of its needed: %d\n', its);
%If exit is equal to zero then the number of its (15) wasnot reached therefore the method converges
if fin == 0
   disp('The point (3.02,-3.03) converges')
end
fprintf('\n\n')

fprintf('QUESTIONS TO BE ANSWERED: \n');

fprintf('From your observations about these 3 points, is the process well-conditioned as a whole?');

fprintf('After analyzing the three points that were given to our team we have concluded this process is ill-conditioned \n')
fprintf('since the condition numbers at the last iteration of P2 and P3 are not approximately 1 (even though P1\n')
fprintf('is close). \n')

fprintf('\n\n')
%(1)What happens if the point’s orbit does, at some step, leave the square [−6, 6] × [−6, 6]?');

fprintf('1. If the point leaves the square, then the Newton-Raphson method will\n')
fprintf('make said point reenter the orbit and then it will converge to a minimum, maximum or saddle point');

%(2) what happens near the crest of the hill in the middle of the square? (observe the plotted function)
fprintf('\n\n')
fprintf('2.Near the middle of the square (near the point 0,0), the algorithm\n')
fprintf('will tend to go to the point (-0.270845, -0.923039), which is a global maximum\n')
fprintf('in the Himmelblau''s function, so all points near (0,0) will tend to go to that point.Also,\n')
disp('there will be found some saddle points near the point so that the curvature change is possible (see the red areas in the plot)');
%(3) what criteria did you choose for stopping?
fprintf('\n\n')
fprintf('3.In our case, we decided to assume that the minimum error for any value during the\n')
fprintf('algorithm should have exactly 6 decimal digits , that is, a tolerance of 10^(-6)\n')
fprintf('The algorithm will continue if the error of the Euclidean norm of the point is bigger.\n')
fprintf('than the required tolerance. It will stop when the current error is smaller or equal to the tolerance.\n')


