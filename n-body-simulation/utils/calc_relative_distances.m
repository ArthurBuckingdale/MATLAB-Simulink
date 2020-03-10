function [relative_distance,grav_potential]=calc_relative_distances(body_information)
%the purpose of this function is two fold. To calculate the relative
%distance for the integrator and for the potential energy calculation.
%it receives the body_information structure and returns the relative
%distance structure array. I have chosen this format of struct arrays since
%there is a good function to unfold and visualise, plus it is much easier
%to track things. I'm thinking this is also the position where we will
%compute the scalar gravitational potential between all the objects. 
%       Inputs: body_information
%                   .name=string, of the body name
%                   .mass=double, mass of the body
%                   .velocity=3x1 vect double, initial velocity
%       Outputs:relative_distance
%                   .body_im_in=string, name of the body were in
%                   .value=double, value of the distance
%                   .body_im_looking_at=string, body that we are looking at
grav_potential(1,1)=zeros();
newtonsg = 6.647e-11;

for i=1:length(body_information)
    for j=1:length(body_information)
        if j~=i
            relative_distance(i).value(j)=norm(body_information(i).position-body_information(j).position);
            relative_distance(i).body_im_in{j}=body_information(i).name;
            relative_distance(i).body_im_looking_at{j}=body_information(j).name;
            grav_potential(i,j)=((newtonsg.*body_information(i).mass.*body_information(j).mass)./(relative_distance(i).value(j)));
        end
    end
end


end