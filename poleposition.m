% Define constants for the transfer function Gp
A = -0.1976068851;  % Constant A in the numerator
K = 198.11;         % Gain K
B = 0.205112;       % Constant B in the denominator
C = 0.557971;       % Constant C in the denominator
D = 1;              % Constant D in the denominator

% Create the transfer function Gp using the given constants
Gp = tf([A*K K],[B C D]);

% Define desired damping ratio and natural frequency
zeta = 0.4559;      % Damping ratio
Wn = 2.3532;        % Natural frequency
P = 10*zeta*Wn;     % Desired pole location

% Symbolically define the PID gains (Kp, Ki, Kd) to solve for
syms Kp Ki Kd

% Define the first equation based on the desired closed-loop pole
eq1 = C/(B - A*K*Kd) + (K*Kd)/(B - A*K*Kd) - (A*K*Kp)/(B - A*K*Kd) == P + 2*zeta*Wn;

% Define the second equation based on the desired closed-loop pole
eq2 = D/(B - A*K*Kd) + (K*Kp)/(B - A*K*Kd) - (A*K*Kp)/(B - A*K*Kd) == 2*Wn*P + Wn^2;

% Define the third equation based on the desired closed-loop pole
eq3 = (A*K*Ki)/(B - A*K*Kd) == P*(Wn^2);

% Solve the system of equations for Kp, Ki, and Kd
S = solve(eq1,eq2,eq3);

% Convert the symbolic solutions to numerical values
X = double(vpa(S.Kp));  % Proportional gain
Y = double(vpa(S.Ki));  % Integral gain
Z = double(vpa(S.Kd));  % Derivative gain

% Create the PID controller transfer function Gc using the calculated gains
Gc = tf([Z X Y],[1 0]);

% Alternative PID controller transfer function for different gains
Gc1 = tf([0.005 0.002 0.02],[1 0]);

Gc2 = tf([0.0035 0.01 0.02],[1 0]);

Gc3 = tf([1 5],[1 10]);

% Plot the step response of the system without feedback
figure(1);
step(Gp);  % step response of open-loop system

% Calculate the closed-loop transfer function with feedback
GK = feedback(Gp,Gc);

GK1 = feedback(Gp,Gc1);

GK2 = feedback(Gp,Gc2);

GK3 = feedback(Gp,Gc3);

% Define a new figure for another plot
figure(2);

% Plot the step response of the new closed-loop systems
step(GK);
hold on;   % hold plot for adding more plots
step(GK1);
step(GK2);
step(GK3);
