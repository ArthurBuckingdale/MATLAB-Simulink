%the purpose of this function is to drive the n-body solveing model for the
%kepler-11 system. This is here to provide a second example, and a more
%interesting one since the planets orbit much faster than in our solar
%system. As well since this was part of the thesis I completed, I've seen
%the results from the c++ code (which is found under the test folder in
%this directory). I have also completed the calculations for all the
%initial conditions. 



%% enter in the mass of each body in the system. 
%for this, I believe they are given as a function of jupiter masses, which
%is why they are all multiplied by it, save for the star ofcourse 
mass(7)=[0.961*1.989e30,0.006*1.898e27,0.009*1.898e27,...
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

%%


    x(0) = -5000000.0; % x1   star
    x(1) = 0.0; % vx1
    x(2) = 0.0 ;% y1
    x(3) = -3.3 ; vy1
    x(4) = -1.01 ;% z1
    x(5) = 1.01 ;% vz1

    x(6) = 1.361000e10*(0.980111);% x2  planet 1   kepler-11b
    x(7) =  9.6554e4*(0.19844);% vx2
    x(8) = 1.361000e10*(-0.19844);% y2
    x(9) =  9.6554e4*(0.980111);% vy2
    x(10) = -6.0 ;% z2
    x(11) = 1.000;% vz2

    x(12) =  1.5993255e10;% x3  planet 2 kepler-11 c
    x(13) = 0.0;% vx3
    x(14) =  0.0;% y3
    x(15) =  8.9301e4;% vy3
    x(16) = 0.000 ;% z3
    x(17) = 0.000 ;% vz3

    x(18)= (-0.95687)*2.315429e10;% x4  planet 3 kepler-11 d
    x(19)= 7.4219e4*(0.290512); % vx4
    x(20)= (-0.290512)*2.315429e10; % y4
    x(21)= 7.4219e4*(-0.95687); % vy4
    x(22)= 0.0;% z4
    x(23)= 0.0; % vz4

    x(24)= 2.91187e10; %x5
    x(25)= (0.0025889)*6.6184e4; %vx5 planet 4 kepler-11 e
    x(26)= -(0.0025889)*2.91187e10; %y5
    x(27)= 6.6184e4;%vy5
    x(28)= 1.0;% z5
    x(29)= 1.0;  %vz5

    x(30)=3.7399000e10*(0.6440731);  %planet 5 kepler-11 f
    x(31)=(-0.764963)*5.8252000e4;
    x(32)=3.7399000e10*(0.7649639);
    x(33)=5.8252000e4*(0.6440731);
    x(34)=1.0;
    x(35)=1.0;

    x(36)=6.97126000e10*(-0.060356);% planet 6 kepler-11 g
    x(37)=4.2824000e4*(0.998176);
    x(38)=6.97126000e10*(-0.998176);
    x(39)=-4.2824000e4*(0.0660356);
    x(40)=1.0;
    x(41)=1.0;
