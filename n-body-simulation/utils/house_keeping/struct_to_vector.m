function x=struct_to_vector(body_information)
%this function is going to shape our data in the format that the ode
%integrator expects it.(remember that this sucks for us because of all the
%different coponents we need to compute). It's going to be much easier to
%compute the values we need when stored in the structure array format(it's
%less than 10 lines of code to compute and is very readable). We can use
%these on more than a few occaisions to help us visualise the data as well.

x = reshape(cat(1,[body_information(:).position; body_information(:).velocity]),[],1);

end