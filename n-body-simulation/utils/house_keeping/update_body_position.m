function body_information=update_body_position(x,body_information)
%the purpose of this function is to take the vector and update the
%information inside of the body_information structure. Remember, the
%integrator is going to inch us along the slope, then we need to
%recalculate the positions to update the accelerations.
j=1;
for i=1:6:((6*length(body_information)))
    body_information(j).position=x(i:i+2);
    body_information(j).velocity=x(i+3:i+5);
    j=j+1;
end