function kinetic_energy=calc_kinetic_energy(body_information)
%the purpose of this function is to compute the kinetic energy of the
%system. This will be responsible for part of the error monitoring for the
%algorithm. Potential energy will be in conjuction with this.
%
%       Inputs: body_information
%                   .name=string, of the body name
%                   .mass=double, mass of the body
%                   .initial_velocity=3x1 vect double, initial velocity
%       Outputs:kinetic_energy
%                   .name=string, name of the body
%                   .value=double, value of the kinetic energy
for i=1:length(body_information)
    kinetic_energy(i).name=body_information(i).name;
    velocity=norm(body_information(i).velocity);
    kinetic_energy(i).value=0.5.*body_information(i).mass.*(velocity.^2);
end
end