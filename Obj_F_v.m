%% Obj_F Object Function_v with vector

%**************************************************************************
% Caculate the Object Function K                                     Mr.Guo
% Ki = |S21(omigaS21)|*|S21(omigaS21)| or Sigma|S11(omigaS11)|*|S11(omigaS11)|
% Input:
% S21        -- Locations of the S21s (lowpass)
% S11        -- Locations of the S11s (lowpass)
% x         -- Current Location (column vector)
% M         -- Coupling Matrix
% Opti_Var  -- Index of the Opti_Coupling Coefficent
% Output:
% Ki        -- Ki(x) (vector)
%**************************************************************************

function Ki = Obj_F_v(x, S21, S11, S22, Freq, delta, M, Opti_Var, W, CF, FBW)

N = length(M)-2;
Ns = length(S11);
Ki = zeros(length(S21)+N+2,1);

G = zeros(N+2,N+2); %-------Port admittance matrix 
G(1,1) = 1; G(N+2,N+2) = 1;
C = eye(N+2,N+2);  %-------Resonantor characteristic admittance matrix 
C(1,1) = 0; C(N+2,N+2) = 0;
L = [x(end-3),x(end-2)];
Theta = [x(end-1),x(end)];
for index = 1:length(x)-4
    M(Opti_Var{index}(1)+1,Opti_Var{index}(2)+1) = x(index);
    M(Opti_Var{index}(2)+1,Opti_Var{index}(1)+1) = M(Opti_Var{index}(1)+1,Opti_Var{index}(2)+1);
end

for index = 1:length(S21)
    invA = inv(G+(Freq(index)+delta)*C+1i*M);
    S21_M = 2*invA(N+2,1);
    TransX = CF*0.5*(FBW*(-1i)*Freq(index)+sqrt(4-Freq(index)*Freq(index)*FBW*FBW));
    Fai21 = L(1)*TransX + L(2)*TransX + Theta(1) + Theta(2);
    Exp21 = exp(1i*Fai21);
    Ki(index) = sqrt(W(index))*abs(Exp21*S21_M-S21(index));

    S11_M = -1+2*invA(1,1);
    Fai11 = L(1)*TransX + L(1)*TransX + Theta(1) + Theta(1);
    Exp11 = exp(1i*Fai11);
    Ki(index+Ns) = sqrt(W(index+Ns))*abs(Exp11*S11_M-S11(index));

    S22_M = -1+2*invA(N+2,N+2);
    Fai22 = L(2)*TransX + L(2)*TransX + Theta(2) + Theta(2);
    Exp22 = exp(1i*Fai22);
    Ki(index+Ns*2) = sqrt(W(index+Ns*2))*abs(Exp22*S22_M-S22(index));
end

end
