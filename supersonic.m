
% This MATLAB program calculates the steady state solution for
%  the Quasi 1D nozzle using Maccormack technique.

% The algorithm for this program is taken from the
% "Computational Fluid Dynamics" book by John D Anderson.
% All the assumed quantities are taken from that book.
%
% We aslo followed the same quantities so that the output values from
% this program can be cross checked with the author's machintosh computer's values.
%
% Note :: The grid points are set to be 30, it is verified that the values are
%  grid independent.
%
% Start of the program

clear all;

% x is the length of the nozzle and each point.
x =0:0.1:3;
%points are the grid points along the length of the nozzle
points = 31;
dx = x(2)-x(1);
%value at the throat
th =16;


% initial grid point values at t = 0
r = 1-0.3146*x;                 %density r
t = 1-0.2314*x;                 %temperature t
v = (0.1+1.09*x).*(sqrt(t));    % velocity v



%assumed quantities
a = 1+2.2*(x-1.5).^2;                   %area equation
C= 0.5;                                 %CFL number
gamma = 1.4;                            % constant ratio
iter = 2000;                            %no of iterations
total_time= 0;                          % total time taken for achieving steady state.
%time = min(C*(dx./(sqrt(t)+v)));
%total_time = total_time+time;
tic
% Iterations loop
for i = 1:iter

% copying the initial data
init_r = r;
init_v = v;
init_t = t;

% calculation of the time
time = min(C*(dx./(sqrt(t)+v)));
total_time = total_time+time;

%predictor step
for j = 2:points-1

% foreward difference terms
vfd = (v(j+1)-v(j))/dx;
rfd = (r(j+1)-r(j))/dx;
tfd = (t(j+1)-t(j))/dx;
afd = (log(a(j+1))-log(a(j)))/dx;


%continuity equation
dr_p(j) = -r(j)*vfd - r(j)*v(j)*afd - v(j)*rfd ;

%momentum equation
dv_p(j) = -v(j)*vfd - (1/gamma)*(tfd +(t(j)/r(j))*rfd ) ;

%energy equation
dt_p(j) = -v(j)*tfd -(gamma-1)*t(j)*(vfd + v(j)*afd) ;

%solution update-predicted values
r(j) = r(j) + dr_p(j)*time;
v(j) = v(j) + dv_p(j)*time;
t(j) = t(j) + dt_p(j)*time;

end
%end of predictor step


%corrector step
for k = 2:points-1

% rearward difference terms
vrd = (v(k)-v(k-1))/dx;
rrd = (r(k)-r(k-1))/dx;
trd = (t(k)-t(k-1))/dx;
ard = (log(a(k))-log(a(k-1)))/dx;

%continuity equation
dr_c(k) = -r(k)*vrd -r(k)*v(k)*ard -v(k)*rrd;

% momentum equation
dv_c(k) = -v(k)*vrd -(1/gamma)*(trd +(t(k)/r(k))*rrd);

%energy equation
dt_c(k) = -v(k)*trd - (gamma-1)*t(k)*(vrd + v(k)*ard) ;

end

% average time step
dr_avg = 0.5*(dr_p + dr_c);
dv_avg = 0.5*(dv_p + dv_c);
dt_avg = 0.5*(dt_p + dt_c);

for l = 2:points-1

%solution update
r(l) = init_r(l) + dr_avg(l)*time;
v(l) = init_v(l) + dv_avg(l)*time;
t(l) = init_t(l) + dt_avg(l)*time;

end

%boundary values - inflow
v(1) = 2*v(2)-v(1) ;

%boundary values - outflow
v(points) = 2*v(points-1) - v(points-2);
r(points) = 2*r(points-1) - r(points-2);
t(points) = 2*t(points-1) - t(points-2);


%non dimensional  parameters
pressure = r.*t;
mach = v./sqrt(t);
mflow = r.*a.*v;


%residuals- storing at throat
r_res(i) = dr_avg(th);
v_res(i) = dv_avg(th);
t_res(i) = dt_avg(th);

%parameters at throat
p_th(i) = pressure(th);
mach_th(i) = mach(th);
r_th(i) = r(th);
t_th(i) = t(th);

%ploting of mass flow rate at throat for  different iterations

figure(1)
hold on

if i ==1
plot(x,mflow,'r','LineWidth',2)
hold on

elseif i ==50
plot(x,mflow,'b','LineWidth',2)
hold on

elseif i ==100
plot(x,mflow,'g','LineWidth',2)
hold on

elseif i == 300
plot(x,mflow,'y','LineWidth',2)
hold on

elseif i == 700
plot(x,mflow,'c','LineWidth',2)
hold on

elseif i ==1400
plot(x,mflow,'LineWidth',2)
hold on

end

%plotting the massflow at different iterations
xlabel("Non Dimensional length");
ylabel("Non Dimensional mass flow rate");
title("Non Dimensional mass flow at different iterations",'FontWeight','bold','fontsize',12);
%legend('iter =1', 'iter = 50', 'iter = 100', 'iter = 300', 'iter = 700', 'iter=1400')



end
legend('1t', '50t', '100t', '300t', '700t', '1400t')
saveas(gcf,'iterations.png');
last = toc ;
%end of iterations loop

%updating the data table after 1400 iterations
table(:,1) = 1:31;
table(:,2) = x;                  % points along the length of the nozzle
table(:,3) = a;                  % area along the each point
table(:,4) = r;                  % dimensionless density values at each grid point
table(:,5) = v;                  % dimensionless velocity values at each grid point
table(:,6) = t;                  % dimensionless temperature values at each grid point
table(:,7) = pressure;           % dimensionless pressure
table(:,8) = mach;               % Mach number
table(:,9) = mflow;              % mass flow rate


total_time= total_time*(2.94*10^-3);

disp("the total time taken for the flow to achive steady state is "+total_time +" s");
disp("total time of computation is "+last +" s");

% ploting of residuals at the throat
figure(2)
plot(r_res)
hold on
plot(v_res)
hold on
plot(t_res)
hold on
xlabel("No of iterations(time steps)")
ylabel(" residuals");
title("Scaled Residuals")
legend('Density','Velocity','Temperature')
grid minor
saveas(gcf,'residuals.png');



%plotting of parameters at throat for different iterations
figure(5)
plot(p_th)
hold on
plot(mach_th)
hold on
plot(r_th)
hold on
plot(t_th)
xlabel('iterations')
ylabel('parameters')
legend('Pressure','Mach Number','Density','Temperature')
title('Variation in parameters at Throat.')
grid on
saveas(gcf,'thoat_parameters.png');

%plotting of parameters along the length of the  nozzle
%mach number
figure(3)
subplot(4,1,3)
plot(x,mach,'b');
xlabel("length");
ylabel("Mach Number");
grid minor
hold on

%pressure ratio
subplot(4,1,4)
plot(x,pressure,'b')
hold on
xlabel('length')
ylabel('Pressure')
grid minor

%density
subplot(4,1,2)
plot(x,r,'b')
hold on
xlabel('length')
ylabel('Density')
grid minor

%temperature
subplot(4,1,1)
plot(x,t,'b')
hold on
xlabel('length')
ylabel('Temperature')
grid minor
title("Qualitive aspects of Quasi - 1D flow");
saveas(gcf,'Qualitative_aspects.png');

%writing data to a xlsx file

writematrix(table,'nozzle_data,xlsx');

% the end of MATLAB program.
