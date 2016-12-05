function dydt = RevKinRxn(t,y,kF,kR)
%RevKinRxn A+B<=>C kinetically
%   kF and kR are rate constants
dydt = zeros(3,1);
dydt(1) = -kF.*y(1).*y(2) + kR.*y(3);
dydt(2) = -kF.*y(1).*y(2) + kR.*y(3);
dydt(3) =  kF.*y(1).*y(2) - kR.*y(3);
end

