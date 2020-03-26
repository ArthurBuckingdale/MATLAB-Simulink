%the purpose of this function is to drive the n-body solveing model for the
%kepler-11 system. This is here to provide a second example, and a more
%interesting one since the planets orbit much faster than in our solar
%system. As well since this was part of the thesis I completed, I've seen
%the results from the c++ code (which is found under the test folder in
%this directory). I have also completed the calculations for all the
%initial conditions. 


clear all
%% enter in the mass of each body in the system. 
%for this, I believe they are given as a function of jupiter masses, which
%is why they are all multiplied by it, save for the star ofcourse 
mass=[0.961*1.989e30,0.006*1.898e27,0.009*1.898e27,...
    0.023*1.898e27,0.025*1.898e27,0.006*1.898e27,0.079*1.898e27];


%% enter in the names for each body
%same as our solar system convention-wise. For extrasolar systems, the star
%is labelled as Kepler-11a, the first planet will be Kepler-11b and so on. 
% not as elegant as names, but oh well. 

names{1}='Kepler-11a';
names{2}='Kepler-11b';
names{3}='Kepler-11c';
names{4}='Kepler-11d';
names{5}='Kepler-11e';
names{6}='Kepler-11f';
names{7}='Kepler-11g';

%% state vectors
% this data could be entered a bit cleaner, but I grabbed it from c++
% code and don't want to rearrange it all. Remember we're dealing with vector
%that go down a colum, not across a row. 


    x(1,1) = -0.0; % x1   star
    x(4,1) = 0.0; % vx1
    x(2,1) = 0.0 ;% y1
    x(5,1) = -3.3 ; %vy1
    x(3,1) = -1.01 ;% z1
    x(6,1) = 1.01 ;% vz1

    x(1,2) = 1.361000e10*(0.980111);% x2  planet 1   kepler-11b
    x(4,2) =  9.6554e4*(0.19844);% vx2
    x(2,2) = 1.361000e10*(-0.19844);% y2
    x(5,2) =  9.6554e4*(0.980111);% vy2
    x(3,2) = -6.0 ;% z2
    x(6,2) = 1.000;% vz2

    x(1,3) =  1.5993255e10;% x3  planet 2 kepler-11 c
    x(4,3) = 0.0;% vx3
    x(2,3) =  0.0;% y3
    x(5,3) =  8.9301e4;% vy3
    x(3,3) = 0.000 ;% z3
    x(6,3) = 0.000 ;% vz3

    x(1,4)= (-0.95687)*2.315429e10;% x4  planet 3 kepler-11 d
    x(4,4)= 7.4219e4*(0.290512); % vx4
    x(2,4)= (-0.290512)*2.315429e10; % y4
    x(5,4)= 7.4219e4*(-0.95687); % vy4
    x(3,4)= 0.0;% z4
    x(6,4)= 0.0; % vz4

    x(1,5)= 2.91187e10; %x5
    x(4,5)= (0.0025889)*6.6184e4; %vx5 planet 4 kepler-11 e
    x(2,5)= -(0.0025889)*2.91187e10; %y5
    x(5,5)= 6.6184e4;%vy5
    x(3,5)= 1.0;% z5
    x(6,5)= 1.0;  %vz5

    x(1,6)=3.7399000e10*(0.6440731);  %planet 5 kepler-11 f
    x(4,6)=(-0.764963)*5.8252000e4;
    x(2,6)=3.7399000e10*(0.7649639);
    x(5,6)=5.8252000e4*(0.6440731);
    x(3,6)=1.0;
    x(6,6)=1.0;

    x(1,7)=6.97126000e10*(-0.060356);% planet 6 kepler-11 g
    x(4,7)=4.2824000e4*(0.998176);
    x(2,7)=6.97126000e10*(-0.998176);
    x(5,7)=-4.2824000e4*(0.0660356);
    x(3,7)=1.0;
    x(6,7)=1.0;
    
%% the orbital periods
%this is quite interesting for the Kepler-11 system. This is due to the
%extremely rapid orbit time compared to human standards. The inner most
%planet orbits the star in 10 days. See the:https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative
%for the location which I retrieved this data. I may in the future create
%an interface to bring this data into the program automatically if there is
%some sort of API for the exoplanet archive. 
orbital_period=[890265.6,1126042.56,1960165.44,2764368.0,4033635.84,10227911.04];

%% subtraction for transit time based on IC's
%since we choose a transit to be a crossing of the x-axis, we're required
%to not only subtract the orbital period, but also add in the initial 
%fraction of the period they've completed already based on their initial
%conditions. Fortunately, we have a way of doing this relatively easily due
%to the way that we've chosen our transits. Since y=0 when we transit, we can 
%easily compute the fraction of an orbit that has passed. I'm almost 100%
%sure there is a better way to perform the below algorithm. I did it late
%at night and it's going to work for now(I can probably query these from the 
%exoplanet archive(maybe not since we've chosen a zero for one))....

init_orbit_fractions=compute_initial_orbit_percentage(x,orbital_period);
disp(init_orbit_fractions)

    
%% choose the debugging boolean
%the debugging option will display all of the calculated outputs from
%various functions. If something is not making sense enable this option to
%see how the data is being passed into the function and used. Check that
%the values make sense. This script is partially a sanity check. We know
%our solar system well and can use it for a validation. 
DEBUG = 1;

%% options
%there will be various options to configure the output of the function or
%tune it to the particular use case that we are having. (no options yet
%code is still being developped.One is the time span. This will yield the
%duration which we integrate the system for. 
plot_xy = 0; %makes a overhead plot of the system
plot_trans_time = 0;
method=2; %1 = ode45, 2=BS

%% set the timespan for integration 
tspan = [0.001 50000];

%% executing the function
%this section will call and execute the fucntion
[transit_timings,integration_error]=newtonian_gravity_nbody(names,mass,x,DEBUG,...
    tspan,plot_xy,plot_trans_time,orbital_period,init_orbit_fractions,method);

disp(init_orbit_fractions)

%% various post processing and such things
%show all the data that comes from this function...
