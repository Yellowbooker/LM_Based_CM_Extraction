%% Jacobian

%**************************************************************************
% Caculate the Jacobian matrix dKi/dMpq                               By ---
% K = Sigma Ki
%   = Sigma|S21(omigaS21)|*|S21(omigaS21)| + Sigma|S11(omigaS11)|*|S11(omigaS11)|
% Input:
% S21      -- Locations of the S21s (lowpass)
% S11      -- Locations of the S11s (lowpass)
% Opti_Var-- [pq]
% M       -- Current Location (matrix)
% G       -- Port Matrix
% C       -- Cap Matrix
% Output:
% Jx      -- dKi/dMpq|x=x0 (matrix)
%**************************************************************************
% clear all

function Jx = Jacobian(S21, S11, S22, Freq, Opti_Var, L, Theta, delta, M, G, C, W, CF, FBW)

N = length(M)-2;
Ns = length(S21);
Jx = zeros(Ns*3,length(Opti_Var)+4);

for index = 1:length(Opti_Var)+4 %    

        for index1 = 1:Ns
            invA = inv(G+(Freq(index1)+delta)*C+1i*M);
            S21_M = 2*invA(N+2,1);
            TransX = CF*0.5*(FBW*(-1i)*Freq(index1)+sqrt(4-Freq(index1)*Freq(index1)*FBW*FBW));
            Fai21 = L(1)*TransX + L(2)*TransX + Theta(1) + Theta(2);
            Exp21 = exp(1i*Fai21);
            if index <= length(Opti_Var)
                p = Opti_Var{index}(1)+1;
                q = Opti_Var{index}(2)+1;
                Jx(index1,index) = sqrt(W(index1))*real(-2i*abs(Exp21*S21_M-S21(index1))/...
                    (Exp21*S21_M-S21(index1))*Exp21*(invA(N+2,p)*invA(q,1)+invA(N+2,q)*invA(p,1)));
            end
            if (index > length(Opti_Var)) && (index <= length(Opti_Var)+2)
                Jx(index1,index) = sqrt(W(index1))*real(1i*abs(Exp21*S21_M-S21(index1))/...
                    (Exp21*S21_M-S21(index1))*Exp21*S21_M*TransX);
            end
            if index > length(Opti_Var)+2
                Jx(index1,index) = sqrt(W(index1))*real(1i*abs(Exp21*S21_M-S21(index1))/...
                    (Exp21*S21_M-S21(index1))*Exp21*S21_M);
            end
            
            S11_M = -1+2*invA(1,1);
            Fai11 = L(1)*TransX + L(1)*TransX + Theta(1) + Theta(1);
            Exp11 = exp(1i*Fai11);
            if index <= length(Opti_Var)
                p = Opti_Var{index}(1)+1;
                q = Opti_Var{index}(2)+1;
                Jx(index1+Ns,index) = sqrt(W(index1+Ns))*real(-4i*abs(Exp11*S11_M-S11(index1))/...
                    (Exp11*S11_M-S11(index1))*Exp11*invA(1,p)*invA(q,1));
            end
            if index == length(Opti_Var)+1
                Jx(index1+Ns,index) = sqrt(W(index1+Ns))*real(2i*abs(Exp11*S11_M-S11(index1))/...
                    (Exp11*S11_M-S11(index1))*Exp11*S11_M*TransX);
            end
            if index == length(Opti_Var)+3
                Jx(index1+Ns,index) = sqrt(W(index1+Ns))*real(2i*abs(Exp11*S11_M-S11(index1))/...
                    (Exp11*S11_M-S11(index1))*Exp11*S11_M);
            end
            
            S22_M = -1+2*invA(N+2,N+2);
            Fai22 = L(2)*TransX + L(2)*TransX + Theta(2) + Theta(2);
            Exp22 = exp(1i*Fai22);
            if index <= length(Opti_Var)
                p = Opti_Var{index}(1)+1;
                q = Opti_Var{index}(2)+1;               
                Jx(index1+Ns*2,index) = sqrt(W(index1+2*Ns))*real(-4i*abs(Exp22*S22_M-S22(index1))/...
                    (Exp22*S22_M-S22(index1))*Exp22*invA(N+2,p)*invA(q,N+2));
            end
            if index == length(Opti_Var)+2
                Jx(index1+Ns*2,index) = sqrt(W(index1+2*Ns))*real(2i*abs(Exp22*S22_M-S22(index1))/...
                    (Exp22*S22_M-S22(index1))*Exp22*S22_M*TransX);
            end
            if index == length(Opti_Var)+4
                Jx(index1+Ns*2,index) = sqrt(W(index1+Ns*2))*real(2i*abs(Exp22*S22_M-S22(index1))/...
                    (Exp22*S22_M-S22(index1))*Exp22*S22_M);
            end
            
            
        end
        

end
end

