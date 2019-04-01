function rA = CorrAB(A,B)
[R,C] = size(A);
[R2,C2] = size(B);
Rw2 = (R2-1)/2;
Cw2 = (C2-1)/2;
Ro = Rw2+1;
Co = Cw2+1;
meanA = zeros(R-2*Rw2,C-2*Cw2);
for i=-Rw2:Rw2
    for j=-Cw2:Cw2
        meanA = meanA + A(1+Rw2+i:R-Rw2+i,1+Cw2+j:C-Cw2+j);
    end
end
meanA = meanA/(R2*C2);
Bhat = B - mean2(B);
SA = zeros(R-2*Rw2,C-2*Cw2);
MA = zeros(R-2*Rw2,C-2*Cw2);
for i=-Rw2:Rw2
    for j=-Cw2:Cw2
        SA = SA + (A(1+Rw2+i:R-Rw2+i,1+Cw2+j:C-Cw2+j)-meanA).*(Bhat(Ro+i,Co+j));
        MA = MA + (A(1+Rw2+i:R-Rw2+i,1+Cw2+j:C-Cw2+j)-meanA).^2;        
    end
end
MA = (MA.*(sum(sum(Bhat.^2)))).^0.5;
rA = SA./MA;

end
