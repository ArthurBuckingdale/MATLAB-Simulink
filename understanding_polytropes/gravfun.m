function gravPotEn = gravfun(n,mass,radius)
%the purpose of this function is to compute the gravitational potential
%enegy that is contained in the polytropic star model. 
G=6.67384*(10.^(-11));
grav_const = -(3./(5-n));
grav_value = G.*((mass.^2)./(radius));
gravPotEn = grav_const.*grav_value;