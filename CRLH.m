clear;
clc;
delta_p = 4; % 4mm < (landa/10)
Er = 4.6;
t = 1.6*10^(-3) ; % the thickness of the substrate = 1.6mm
tan_alpha = 0.02;
fx = 2*10^(9):10^(8):2.8*10^(9);
f_point = 2.4*10^9; % Design Frequency = 2.4GHz
S21x = [];
S21x;
for i = fx
    f0 = i;
    c = 30000000; % speed of light = 300000000 m/s
    omega = 2*3.14*f0; % omega
    omega_point = 2*3.14*f_point; % omega
    %%% Desired phase shift (Degree) = +10 
    %%% Computed Value From (12) and (13) for single unit cell 
    %%% CRLH phase shifter
    Ll = 6.043^(-9); %% unit[nH]
    Cl = 2.417^(-12); %% unit[pF]
    W = 3*10^(-3); % the conductor width = 3mm
    d = 4*10^(-3); % d = 4mm
    Z0 = sqrt(Ll / Cl); % unit[ohm]
    Y0 = 1/Z0;
    k0 = (2*3.14*f0)/c; %The free space wave number
    k0_point = (2*3.14*f_point)/c; 
    Ee = (Er + 1)/2 + ((Er - 1)/2)*(1/sqrt(1 + (12*t/W)));
    betha_TL = sqrt(Ee)*k0;
    betha_TL_point = sqrt(Ee)*k0_point;
    A1 = cos(betha_TL*d);
    A1_point = cos(betha_TL_point*d);
    B1 = 1i*Z0*sin(betha_TL*d);
    B1_point = 1i*Z0*sin(betha_TL_point*d);
    C1 = 1i*Y0*sin(betha_TL*d);
    C1_point = 1i*Y0*sin(betha_TL_point*d);
    D1 = A1;
    D1_point = A1_point;
    Z1 = rdivide(1,(1i*omega*2*Cl));
    Z1_point = rdivide(1,(1i*omega_point*2*Cl));
    Z2 = Z1;
    Z2_point = Z1;
    Z3 = 1i*omega*Ll;
    Z3_point = 1i*omega_point*Ll;
    A2 = 1 + Z1/Z3;
    A2_point = 1 + Z1_point/Z3_point;
    B2 = Z1 + Z2 + Z1*Z2/Z3;
    B2_point = Z1_point + Z2_point + Z1_point*Z2_point/Z3_point;
    C2 = 1/Z3;
    C2_point = 1/Z3_point;
    D2 = 1 + Z2/Z3;
    D2_point = 1 + Z2_point/Z3_point;
    array_1 = [A1 B1;C1 D1];
    array_1_point = [A1_point B1_point;C1_point D1_point];
    array_2 = [A2 B2;C2 D2];
    array_2_point = [A2_point B2_point;C2_point D2_point];
    array = array_1.*array_2;
    array_point = array_1_point.*array_2_point;
    pi_PRH = betha_TL * d;
    pi_PRH_point = betha_TL_point * d;
    pi_PLH = -1/(omega*sqrt(Ll*Cl));
    pi_PLH_point = -1/(omega_point*sqrt(Ll*Cl));
    pi_CRLH = pi_PRH + pi_PLH;
    pi_CRLH_point = pi_PRH_point + pi_PLH_point;
    S21 = 2/(array(1) + array(2)/Z0 + array(3)*Z0 + array(4));
    S21_point = 2/(array_point(1) + array_point(2)/Z0 + array_point(3)*Z0 + array_point(4));
    IL = -20*log10(abs(S21));
    IL_point = -20*log10(abs(S21_point));
    S21x(end+1,1) = S21;
end
plot(fx,S21x);
xlabel('Frequency(Hz)')
ylabel('Phase of S21 (deg)')
title('CRLH Figure')
xtickformat('%.1f')
hold on
plot(f_point,S21_point,'r*')
%%axis([0 120 0 4000])
grid on




