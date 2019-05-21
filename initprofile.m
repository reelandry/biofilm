function UU = initprofile(XX, YY, Li,  C, r, xy)

UU = zeros(size(XX));
for ji = 1:Li
UU = UU + C(ji)*exp(-r(ji)*((XX-xy(ji,1)).^2+(YY-xy(ji,2)).^2));
end

return 



