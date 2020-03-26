function init_orbit_fractions=compute_initial_orbit_percentage(x,orbital_period)
%the purpose of this routine is to compute the fraction of each orbit
%that's been completed when the routine starts. This is so we can correctly
%obtain the transit timings from here on out.

[~,mm]=size(x);
init_orbit_fractions(length(orbital_period))=zeros();
for i=2:mm
    angle=atan(x(2,i)/x(1,i)); %this will provide quadrant and location
    disp(angle)
    if x(2,i)>0 %this will provide front or back half
        if angle > 0
            fraction_of_orbit=((2*pi)-angle)/(2*pi);
        elseif angle<0
            fraction_of_orbit=(pi-angle)/(2*pi);
        else
            fraction_of_orbit=0; 
        end
    else        
        if angle > 0 
            fraction_of_orbit=(pi-angle)/(2*pi);
        elseif angle<0
            fraction_of_orbit=((pi/2)+angle)/(2*pi);
        else
            fraction_of_orbit=0; %on the x axis now
        end
    end
    init_orbit_fractions(i)=fraction_of_orbit;
end
