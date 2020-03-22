function dxdt=distance_derivatives_func(t,x,body_information)
%the purpose of this function is to compute the components that need to be
%passed into the differential equation solver this will be iterated over to
%inch forward along the slope to determine the position as a function of
%time. 
%       Inputs: body_information
%                   .name=string, of the body name
%                   .mass=double, mass of the body
%                   .velocity=3x1 vect double, velocity
%                   .positiont=3x1 vect double, position of the body 
%               accel_of_gravity
%                   .separation_distance: N-1 double vect with seperation
%                   distances from each body
%                   .acceleration: 3x(N-1) double vect, has the
%                   acceleration due to gravity for (x,y,x) of each body
%                   .body_im_in:N-1 cell string vect, contains the name of
%                   the body we're standing in
%                   .body_im_looking_at: N-1 cell string vect, containst
%                   the list of bodies that i'm looking at.
%       Outputs:dxdt(), vector that contains the values of the
%               accelerations and the velocities for each body in the system. 

body_information=update_body_position(x,body_information);
accel_of_gravity=calc_grav_accel(body_information);
j=1;

for i=1:6:(6*length(body_information))
    dxdt(i)= body_information(j).velocity(1);
    dxdt(i+1)=body_information(j).velocity(2);
    dxdt(i+2)=body_information(j).velocity(3);
    dxdt(i+3)=sum(accel_of_gravity(j).acceleration(1,:));
    dxdt(i+4)=sum(accel_of_gravity(j).acceleration(2,:));
    dxdt(i+5)=sum(accel_of_gravity(j).acceleration(3,:));
    j=j+1;
end
dxdt=dxdt';
end


