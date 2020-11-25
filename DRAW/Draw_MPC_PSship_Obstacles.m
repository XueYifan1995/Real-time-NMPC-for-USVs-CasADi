function Draw_MPC_PSship_Obstacles (t,xx,xx1,u_cl,xs,N,rob_diam,obs_x,obs_y,obs_diam)

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;

%--------------------------------------------------------------------------
%-----------------------Simulate robots -----------------------------------
%--------------------------------------------------------------------------
x_r_1 = [];
y_r_1 = [];



r = rob_diam/2;  % obstacle radius
ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

r = obs_diam/2;  % obstacle radius
xp_obs=r*cos(ang);
yp_obs=r*sin(ang);

figure(500)
% Animate the robot motion
%figure;%('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

tic
for k = 1:size(xx,2)
    h_t = 4.88; w_t=2.44; % triangle parameters
    
    x1 = xs(4); y1 = xs(5); th1 = xs(6);
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot reference state
    hold on;
    x1 = xx(4,k,1); y1 = xx(5,k,1); th1 = xx(6,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
    if k < size(xx,2) % plot prediction
        plot(xx1(1:N,4,k),xx1(1:N,5,k),'r--*')
        for j = 2:N+1
            plot(xx1(j,4,k)+xp,xx1(j,5,k)+yp,'--r'); % plot robot circle
        end
    end
    
    fill(x1_tri, y1_tri, 'r'); % plot robot position
    
    plot(x1+xp,y1+yp,'--r'); % plot robot circle
    
    plot(obs_x+xp_obs,obs_y+yp_obs,'--b'); % plot robot circle    
    
    hold off
    %figure(500)
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    axis([0 70 0 70])
    pause(0.0)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end
toc
% close(gcf)
%viobj = close(aviobj)
% video = VideoWriter('exp.avi','Uncompressed AVI');
% 
%  video = VideoWriter('exp.avi','Motion JPEG AVI');
%  video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
%  open(video)
%  writeVideo(video,F)
%  close (video)

U_pre = xx(1,2:end);V_pre = xx(2,2:end); R_pre = xx(3,2:end);
figure
subplot(511)
plot(t',U_pre,'linewidth',1.5);xlabel('time (s)');ylabel('u (m/s)');grid on
subplot(512)
plot(t',V_pre,'linewidth',1.5);xlabel('time (s)');ylabel('v (m/s)');grid on
subplot(513)
plot(t',R_pre,'linewidth',1.5);xlabel('time (s)');ylabel('r (rad/s)');grid on
subplot(514)
stairs(t,u_cl(:,1),'r','linewidth',1.5); xlabel('time (s)');ylabel('T_p (N)');grid on
subplot(515)
stairs(t,u_cl(:,2),'r','linewidth',1.5); xlabel('time (s)');ylabel('T_s (N)');grid on
