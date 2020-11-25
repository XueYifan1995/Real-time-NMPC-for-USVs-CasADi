% point stabilization + Single shooting
clear all
close all
clc

import casadi.*

T = 0.12; % sampling time [s]
N = 10; % prediction horizon
rob_diam = 4.88; %visualization

Tp_max = 400;  Tp_min = -400;
Ts_max = 400;  Ts_min = -400;

u=SX.sym('u'); v = SX.sym('v'); r = SX.sym('r');
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');
states = [u;v;r;x;y;psi]; n_states = length(states);

Tp = SX.sym('Tp');Ts = SX.sym('Ts');
controls = [Tp;Ts]; n_controls = length(controls);
% mathetical model

rhs = [-1.1391*u+0.0028*(Tp+Ts)+0.6836;
       0.0161*v-0.0052*r+0.002*(Tp-Ts)*2.44/2+0.0068;
       8.2861*v-0.9860*r+0.0307*(Tp-Ts)*2.44/2+1.3276;  
       u*cos(psi)-v*sin(psi);
       u*sin(psi)+v*cos(psi);
       r]; % system r.h.s plus part

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states); %stares and ref
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X',n_states,(N+1));
% A Matrix that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(6,6); Q(1,1)=0.1; Q(2,2)=0;Q(3,3)=0.1; Q(4,4) = 100;Q(5,5) = 100;Q(6,6) = 0.1; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.001;  R(2,2) = 0.002;% weighing matrices (controls)

st = X(:,1); %initial state
g = [g;st-P(1:6)]; %initial condition constrains

for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con; % calculate obj
%     %4-RK
%     k1 = f(st, con);   % new 
%     k2 = f(st + T/2*k1, con); % new
%     k3 = f(st + T/2*k2, con); % new
%     k4 = f(st + T*k3, con); % new
%     st_next=st +T/6*(k1 +2*k2 +2*k3 +k4); % new   
    st_next = X(:,k+1);
    f_value  = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end

% Add constraints for collision avoidance
obs_x = 40; % meters
obs_y = 20; % meters
obs_diam = 4; % meters

for k = 1:N+1   % box constraints due to the map margins
    g = [g ; -sqrt((X(4,k)-obs_x)^2+(X(5,k)-obs_y)^2) + (rob_diam/2 + obs_diam/2)];
end

% make the decision variables one column vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:6*(N+1)) = 0; % equality constraints
args.ubg(1:6*(N+1)) = 0; % equality constraints

args.lbg(6*(N+1)+1 : 6*(N+1)+ (N+1)) = -inf; % inequality constraints
args.ubg(6*(N+1)+1 : 6*(N+1)+ (N+1)) = 0; % inequality constraints

args.lbx(1:6:6*(N+1),1) = -2; %state u lower bound
args.ubx(1:6:6*(N+1),1) = 2.5; %state u upper bound
args.lbx(2:6:6*(N+1),1) = -inf; %state v lower bound
args.ubx(2:6:6*(N+1),1) = inf; %state v upper bound
args.lbx(3:6:6*(N+1),1) = -inf; %state r lower bound
args.ubx(3:6:6*(N+1),1) = inf; %state r upper bound

args.lbx(4:6:6*(N+1),1) = -100; %state x lower bound
args.ubx(4:6:6*(N+1),1) = 100; %state x upper bound
args.lbx(5:6:6*(N+1),1) = -100; %state y lower bound
args.ubx(5:6:6*(N+1),1) = 100; %state y upper bound
args.lbx(6:6:6*(N+1),1) = -inf; %state theta lower bound
args.ubx(6:6:6*(N+1),1) = inf; %state theta upper bound

% input constraints
args.lbx(6*(N+1)+1:2:6*(N+1)+2*N,1) = Ts_min; %v lower bound
args.ubx(6*(N+1)+1:2:6*(N+1)+2*N,1) = Tp_max; %v upper bound
args.lbx(6*(N+1)+2:2:6*(N+1)+2*N,1) = Ts_min; %omega lower bound
args.ubx(6*(N+1)+2:2:6*(N+1)+2*N,1) = Tp_max;  %omega upper bound

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0;0;0; 0 ;    0 ; 0.0];    % initial condition.
xs = [2.5;0;0;  75 ; 30; 20*pi/180]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);  % 2 control inputs 
X0 = repmat(x0,1,N+1)';  % initialization of the states decision variables

sim_tim = 24; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),1.5) > 1e-2 && mpciter < sim_tim / T)
% while(norm((x0-xs),2) > 2)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0  = [reshape(X0',6*(N+1),1);reshape(u0',2*N,1)];
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p); 
    %toc
    u = reshape(full(sol.x(6*(N+1)+1:end))',2,N)'; % get controls only from the solution
    xx1(:,1:6,mpciter+1)= reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+2) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift_control2ship(T, t0, x0, u,f); % get the initialization of the next optimization step
    xx(:,mpciter+2) = x0;  
    X0 = reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
    mpciter
    mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop);
ss_error = norm((x0-xs),2)
average_mpc_time = main_loop_time/(mpciter+1)

U_pre = xx(1,2:end);V_pre = xx(2,2:end); R_pre = xx(3,2:end);
X_pre = xx(4,:);Y_pre = xx(5,:);

Draw_MPC_PSship_Obstacles (t,xx,xx1,u_cl,xs,N,rob_diam,obs_x,obs_y,obs_diam)