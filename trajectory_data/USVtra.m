%https://www.researchgate.net/publication/271549212_Stabilizing_NMPC_of_wheeled_mobile_robots_using_open-source_real-time_software
dt = 0.1;
S = 4;
tol = 60;
t = linspace(0,tol,tol/dt);

x_r1 = (2*sin(0.25*t))*S-1;
y_r1 = (2*cos(0.25*t))*S-6.2;
psi_r1 = 0.*t;
u_r1 = (x_r1(1:tol/dt-1)-x_r1(2:tol/dt))/dt;
v_r1 = (y_r1(1:tol/dt-1)-y_r1(2:tol/dt))/dt;
r_r1 =  0.*t;
figure(1)
plot(x_r1,y_r1);
ref_usv_circle =  [[0,u_r1];[0,v_r1];r_r1;x_r1;y_r1;psi_r1;];

%1.8 2
x_r2 = 1.7*sin(0.3*t)*S-0.5;
y_r2 = 2.4*cos(0.15*t)*S-6;
psi_r2 = 0.*t;
u_r2 = (x_r2(1:tol/dt-1)-x_r2(2:tol/dt))/dt;
v_r2 = (y_r2(1:tol/dt-1)-y_r2(2:tol/dt))/dt;
r_r2 =  0.*t;
figure(2)
plot(x_r2,y_r2);axis([-20 20 -25 15])
ref_usv_8 = [[0,u_r2];[0,v_r2];r_r2;x_r2;y_r2;psi_r2;];

tol = 60;
dt =0.15;
t = linspace(0,tol,tol/dt);

x_r3 = zeros(1,tol/dt);
y_r3 = 60-2*t;
y_r4 = 2*t;
psi_r3 = zeros(1,tol/dt);
u_r3 = 2* ones(1,tol/dt);
v_r3 = zeros(1,tol/dt);
r_r3 = zeros(1,tol/dt);
ref_usv_line = [u_r3;v_r3;r_r3;x_r3;y_r3;psi_r3;x_r3;];
ref_usv_line2 = [u_r3;v_r3;r_r3;x_r3;y_r4;psi_r3;x_r3;];

save ref_usv_circle ref_usv_circle
save ref_usv_8 ref_usv_8
save ref_usv_line ref_usv_line
save ref_usv_line2 ref_usv_line2