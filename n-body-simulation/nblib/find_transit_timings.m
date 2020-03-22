function transit_time=find_transit_timings(t1,y1,body_information)
%this section will be measuring the point where we're seeing a transit. for
%the sake of this exercise, a transit will occurr when the body is crossing
%the positive x axis. The integrator is not guarenteed to provide a value
%exactly at the zero point for its orbit. For this reason,we must perform a
%small interpolation to obtain its exact transit timing. We don't want to
%fit a sinusoid curve since any elliptic behaviour will throw off these
%values. The routine will parse through the output from the integrator and
%find locations when y(t)==negative and y(t+1)==positive(since our orbits
%are counterclockwise, they will be passing the positive x axis from the
%bottom) we can change this is deemed necessary. We can also improve this
%by removing the linear integrator.
%       input:  t1:matrix, this is what's returned by the integrator. it
%               contains the time stamps for all the observations in the
%               follwing variable
%               y1:matrix, these are the variables returned from the
%               integration. They are statevectors in the same format that
%               the integrator returns(see the distance_derivative_func)
%               body_information: body_information, struct
%                   .name=string, of the body name
%                   .mass=double, mass of the body
%                   .initial_velocity=3x1 vect double, initial velocity
%       outputs:transit_times: struct:
%                   .value:double vect, contains all measured transit times
%                   .body: string, the body that these times are associated
%                   with



fprintf('the length of the time array is %d, the length of the postion array is %d \n',length(t1),length(y1))
[~,vv]=size(y1);
fprintf('the number of objects in this system is %d \n',vv./6)
j=1;
for ss=1:6:vv
    k=1;
    for uu=1:(length(t1)-1)
        if y1(uu,ss+1)<0 && y1(uu+1,ss+1)>0
            x = [t1(uu) t1(uu+1)];
            y = [y1(uu,ss+1) y1(uu+1,ss+1)];
            lin_fit = polyfit(x,y,1);
            transit_time(j).value(k)=-(lin_fit(2)./lin_fit(1));
            transit_time(j).body=body_information(j).name;
            k=k+1;
        end
    end
    j=j+1;
end