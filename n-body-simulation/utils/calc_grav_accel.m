function accel_of_gravity=calc_grav_accel(body_information)
%the purpose of this function is to copmute the first and second
%derivative that are a result of the force of gravity. It will also
%store the body we're on and the body we're looking at.
%       Inputs: body_information
%                   .name=string, of the body name
%                   .mass=double, mass of the body
%                   .velocity=3x1 vect double, velocity
%                   .positiont=3x1 vect double, position of the body 
%       Outputs:accel_of_gravity
%                   .separation_distance: N-1 double vect with seperation
%                   distances from each body
%                   .acceleration: 3x(N-1) double vect, has the
%                   acceleration due to gravity for (x,y,x) of each body
%                   .body_im_in:N-1 cell string vect, contains the name of
%                   the body we're standing in
%                   .body_im_looking_at: N-1 cell string vect, containst
%                   the list of bodies that i'm looking at.

%define the format before entering the parfor loop
% accel_of_gravity(1).separation_distance(1)=0.0;
%accel_of_gravity(1).acceleration(1:3,2)=zeros();
newtonsg = 6.647e-11;
% accel_of_gravity(1).body_im_in{1}='str';
% accel_of_gravity(1).body_im_looking_at{1}='str';

%entering the loop
for i=1:length(body_information)
    for j=1:length(body_information)
        if j~=i
            separation_distance=(norm(body_information(i).position-body_information(j).position)).^3;       
            accel_of_gravity(i).acceleration(1:3,j)=(newtonsg).*-((body_information(j).mass).*([body_information(i).position]-[body_information(j).position]))./(separation_distance);
            accel_of_gravity(i).body_im_in{j}=body_information(i).name;
            accel_of_gravity(i).body_im_looking_at{j}=body_information(j).name;
        end
    end
end
end