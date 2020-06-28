function dydt = odefcn(t,y,n)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -y(1).^n-(2/t)*y(2);
end