%% Levenberg¨CMarquardt Method (LM)

%**************************************************************************
% LM Search to Minimize the Cost Function                          By ---
% LM (The Levenberg¨CMarquardt Method) is an second-order convergence gradient 
% based method.
%--------------------------------------------------------------------------
% Parameters:
% k          -- Iteration time
% kmax       -- Max iteration time
% u          -- Damping parameter
% v          -- Increase factor of u
% t          -- Initial u factor
% x          -- Independent variables (column vector)
% J          -- Jacobian matrix
% A          -- Approximate Hession J(x)T * J(x) 
% g          -- Gradient J(x)T * K(x)
% hlm        -- LM step hlm = -(A+u*I)T * g
% GR         -- Gain ratio
% N          -- Order of the filter
% Opti_Var   -- Index of the Opti_Coupling Coefficent
% S21         -- Lowpass frequency of the transmission zeros
% S11         -- Lowpass frequency of the reflection zeros
% Res        -- Results of the object function K(x)
% er1        -- End loop flag 1 
% er2        -- End loop flag 2
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Update Time1: 2021.01.25
%**************************************************************************

clear variables
tic
%--------------------------------------------------------------------------
% Initial loop parameter
t = 0.0000001;
kmax = 200;
er1 = 0.0000000001; er2 = 0.0000000001; 
% Initial search location
x = [0 0 0 0 0 1 1 0.865 0.6357 0.6357 0.865 0 0 0 0].'; % Coupling Coeffient
% L = [0,0]; % 2*pi/Vg*L
% Theta = [0,0]; % Phase load
% x = [0.02 0.2 0.1 -0.1 0.3 1.1 1 0.665 0.7357 0.6357 1.365].'; 
% Index of the Opti_Coupling Coefficent
Opti_Var = {[1,1],[2,2],[3,3],[4,4],[5,5],[0,1],[1,2],[2,3],[3,4],[4,5],[5,6]};
Q = 1700;

%--------------------------------------------------------------------------
% Initial Matrix 
% ms1 = 1.024273084; m4l = 1.024273084;
N = 5;
M = zeros(N+2,N+2);
% M(1,2) = ms1; M(N+2,N+1) = m4l;
% M(2,1) = ms1; M(N+1,N+2) = m4l;

[Freq, S11, S21, S22] = Reads2p();
[CF,BW,n,nz] = Parasinput();
NP = length(Freq);

W = ones(NP*3); % weight of f(x)
% for index = 1:NP*3
%     if index <= NP
%         W(index) = 1/sqrt(abs(S21(index)));
%     end
%     if index > NP && index <= NP*2
%         W(index) = 1/sqrt(abs(S11(index-NP)));
%     end
%     if index > NP*2
%         W(index) = 1/sqrt(abs(S22(index-NP*2)));
%     end
% end

for index = 1:NP*3
    if index <= NP
        W(index) = (1/abs(S21(index)))^0.5;
    end
    if index > NP && index <= NP*2
        W(index) = (1/abs(S11(index-NP)))^0.8;
    end
    if index > NP*2
        W(index) = (1/abs(S22(index-NP*2)))^0.8;
    end
end


FBW = BW/CF;
w_low=(Freq./CF-CF./Freq)/FBW;
s_low=1i*w_low;
CF = CF./1000000000;
for index = 1:length(x)-4
    M(Opti_Var{index}(1)+1,Opti_Var{index}(2)+1) = x(index);
    M(Opti_Var{index}(2)+1,Opti_Var{index}(1)+1) = M(Opti_Var{index}(1)+1,Opti_Var{index}(2)+1);
end
delta = 1/Q/FBW;
% % Reflection zeros and transmission zeros
% S11 = [-1.5i 1.5i];
% S11 = [-0.972188i -0.740723i -0.281297i 0.281297i 0.740723i 0.972188i];

G = zeros(N+2,N+2); %-------Port admittance matrix 
G(1,1) = 1; G(N+2,N+2) = 1;
C = eye(N+2,N+2);  %-------Resonantor characteristic admittance matrix 
C(1,1) = 0; C(N+2,N+2) = 0;

%--------------------------------------------------------------------------
% Begin loop
% Initial loop parameters
k = 0;
J = Jacobian(S21, S11, S22, s_low, Opti_Var,[x(end-3),x(end-2)], [x(end-1),x(end)], delta, M, G, C, W, CF, FBW);
A = J.'*J;
g = J.'*Obj_F_v(x, S21, S11, S22, s_low, delta, M, Opti_Var ,W, CF, FBW);
u = t*max(diag(A));
v = 2; 
found = norm(g,inf)-er1;
Res = zeros(kmax,1);

while (found > 0) && (k < kmax)
    k = k+1;
    hlm = -1*(A+u*A.*eye(length(A)))\g;
%     hlm = -1*(A+u*eye(length(A)))\g;
    if norm(hlm) <= er2*(norm(x)+er2)
        found = -1;
    else
        xnew = x + hlm;
        GR = (Obj_F(x, S21, S11, S22, s_low, delta, M, Opti_Var, W, CF, FBW)-...
            Obj_F(xnew, S21, S11, S22, s_low, delta, M, Opti_Var, W, CF, FBW))...
            /(0.5*hlm.'*(u*hlm-g));
        if GR > 0
            x = xnew;
            % Updata Coupling Matrix
            for index1 = 1:length(x)-4
                M(Opti_Var{index1}(1)+1,Opti_Var{index1}(2)+1) = x(index1);
                M(Opti_Var{index1}(2)+1,Opti_Var{index1}(1)+1) =...
                    M(Opti_Var{index1}(1)+1,Opti_Var{index1}(2)+1);
            end
            J = Jacobian(S21, S11, S22, s_low, Opti_Var,[x(end-3),x(end-2)], [x(end-1),x(end)], delta, M, G, C, W, CF, FBW);
            A = J.'*J;
            g = J.'*Obj_F_v(x, S21, S11, S22, s_low, delta, M, Opti_Var, W, CF, FBW);
            found = norm(g,inf)-er1;
            u = u*max([1/3,1-(2*GR-1)^(3)]);
            v = 2;
            Res(k) = Obj_F(x, S21, S11, S22, s_low, delta, M, Opti_Var, W, CF, FBW);
        else
            u = u*v;
            v = 2*v;
            Res(k) = Obj_F(x, S21, S11, S22, s_low, delta, M, Opti_Var, W, CF, FBW);
        end
    end
end
%--------------------------------------------------------------------------
figure(1)
plot(1:k, 10*log10(Res(1:k)),'linewidth',2);
grid on 
%--------------------------------------------------------------------------
% Draw the opti_coupling matrix
F_str = -8; F_sto = 8;
Number = 1500;
R(1,1)=1;
R(n+2,n+2)=1;
W=eye(n+2,n+2);
W(1,1)=0;
W(n+2,n+2)=0;
for k=1:Number
    wp(k)=F_str+(F_sto-F_str)*k/Number;
    TransX = CF*0.5*(FBW*wp(k)+sqrt(4+wp(k)*wp(k)*FBW*FBW));
    Fai21 = x(end-3)*TransX + x(end-2)*TransX + x(end-1) + x(end);
    Exp21 = exp(1i*Fai21);
    Fai11 = x(end-3)*TransX + x(end-3)*TransX + x(end-1) + x(end-1);
    Exp11 = exp(1i*Fai11);
    Fai22 = x(end-2)*TransX + x(end-2)*TransX + x(end) + x(end);
    Exp22 = exp(1i*Fai22);
    A=R+W.*(1i*wp(k)+delta)+(1i).*M;
    AP=inv(A);
    NS21(k)=2*AP(n+2,1)*Exp21;
    NS11(k)=(-1+2*AP(1,1))*Exp11;
    NS22(k)=(-1+2*AP(n+2,n+2))*Exp22;
end

figure(2)
plot(w_low, 20*log10(S11),'r:', w_low, 20*log10(S21),'b:',...
    wp, 20*log10(NS11), wp, 20*log10(NS21), 'LineWidth',2);
legend('S11 3D', 'S21 3D', 'S21 CM','S21 CM','Location', 'NorthWest');
title('S-Parameters Mag');
grid on
figure(3)
plot(w_low,angle(S11),'r:', w_low, angle(S21),'b:',...
    wp, angle(NS11), wp, angle(NS21), 'LineWidth',2);
grid on 
legend('S11 3D', 'S21 3D', 'S21 CM','S21 CM','Location', 'NorthWest');
title('S-Parameters Phase');
figure(4)
plot(w_low, 20*log10(S22),'r:', wp, 20*log10(NS22), 'LineWidth',2);
title('S-Parameters Phase');
grid on
disp([x(end-3),x(end-2)])
disp([x(end-1),x(end)])
%--------------------------------------------------------------------------
toc

