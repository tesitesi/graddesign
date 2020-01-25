% ���P�b�g�G���W���Đ���p�v�Z

gamma_g = 1.125  %�K�X��M��
gamma_c = 1.66      %��p�ܔ�M��
Pr_g = 4 * gamma_g / (9 * gamma_g - 5)    %�K�X�v�����g����
Tcin = 112  %��p�ܓ������x [K]
Cpg = 8000  %�R�ăK�X�舳��M[J/kg K]
Cpc = 6000   %��p�ܒ舳��M [J/kg K]
Rt = 0.05788203 %�X���[�g���a [m]
Dt = Rt * 2     %�X���[�g�a [m]
At = Rt^2 * pi  %�X���[�g�ʐ�
Tc = 3510.96    %�`�����o���x [K]
g = 9.80655 %�d�͉����x[m/s^2]
Pc = 4.0 *10^6  %�`�F���o����[Pa]
Gc = 6.8182 %��p�ܗ���[kg/s]
Gg = 30.0   %���i�ܗ���[kg/s]
c = Pc * At / Gg   %�����r�C���x[m/s]
lambda = 390    %�M�`����[W/mK]
lambda_cond = 0.189  %��p�ܔM��R[W/m*K]
mu_c = 0.000123  %��p�ܔS���W��[Pa*s]
rho_c = 426.32
Pr_c = mu_c * Cpc / lambda_cond    %��p�܃v�����g����

deltaP = 0

x = cccalcS5(:,1)    %���������W
r = cccalcS5(:,2)   %�m�Y�����a
A = cccalcS5(:,3)   %�m�Y���ʐ�
rho = cccalcS5(:,4) %�R�ăK�X���x
T = cccalcS5(:,16)   %�R�ăK�X���x
p = cccalcS5(:,6)   %�R�ăK�X����
M = cccalcS5(:,7)   %�R�ăK�X�}�b�n��
ug = cccalcS5(:,8)  %�R�ăK�X���x
b = cccalcS5(:,9)   %��p�ܗ��H����
t = cccalcS5(:,10)  %��p�ܗ��H�ԕ�
e = cccalcS5(:,11)  %��p�Ǔ���
a = cccalcS5(:,12)  %��p�ܗ��H��
Ac = cccalcS5(:,13) %��p�ܗ��H�ʐ�
uc = cccalcS5(:,14)*426/220 %��p�ܗ���
Dc = cccalcS5(:,15) %��p�ܗ��H�������a
lambda = cccalcS5(:,20)
Dp = cccalcS5(:,21)

throat = 745   %�X���[�g�ʒu
L = size(x)
L = L(1)

Twg = cccalcS5(:,17)
Twg(L)=201
Twc = cccalcS5(:,18)
Tcool = cccalcS5(:,19)
Tcool = vertcat(Tcool,Tcin)

mu = 46.6*10^(-10)*(2.20462*21.25)^0.5*(1.8*T(throat))^0.6
hg_0 = 1635339.84*9/5*0.026/(Dt*39.3701)^0.2*mu^0.2*(gamma_g*1544/(22.2*2.20462)/(gamma_g-1)/778)/Pr_g^0.6*((0.000145038*Pc)*32.2/(c*3.28))^0.8*(Dt/Rt)^0.1

for i = 0:L-1
    deltaT = 1
    exT = 1
    exxT = 1
    while abs(deltaT) > 0.0001
        exxT = exT
        exT = deltaT
        exTwg = Twg(L-i)
        if deltaT > 0
            Twg(L-i) = Twg(L-i) + 0.1
        elseif deltaT < 0
            Twg(L-i) = Twg(L-i) - 0.1
        end
    
        sigma = 1/((0.5*Twg(L-i)/Tc*(1+(gamma_g-1)/2*M(L-i).^2)+0.5)^0.68*(1+(gamma_g-1)/2*M(L-i).^2)^0.12)
        hg = hg_0*(At/A(L-i))^0.9*sigma

        q = hg*(T(L-i) - Twg(L-i))
        deltaT1 = q*0.001*2*pi*r(L-i)/(Gc*Cpc)
        Twc(L-i) = Twg(L-i) - e(L-i)*q/lambda(L-i)

        Re = rho_c*uc(L-i)*Dc(L-i)/mu_c
        hc = lambda_cond/Dc(L-i)*0.027*Re^0.8*Pr_c^(1/3)
        Tcool(L-i) = Twc(L-i)-q/hc

        deltaT2 = Tcool(L-i)-Tcool(L+1-i)
        deltaT = deltaT1 - deltaT2
        
        
        if exxT == deltaT && exT*exxT < 0 && exxT ~= 1
            if abs(exT) < abs(deltaT)
                Twg(L-i) = exTwg
                sigma = 1/((0.5*Twg(L-i)/Tc*(1+(gamma_g-1)/2*M(L-i).^2)+0.5)^0.68*(1+(gamma_g-1)/2*M(L-i).^2)^0.12)
                hg = hg_0*(At/A(L-i))^0.9*sigma

                q = hg*(T(L-i) - Twg(L-i))
                Twc(L-i) = Twg(L-i) - e(L-i)*q/lambda(L-i)

                Re = rho_c*uc(L-i)*Dc(L-i)/mu_c
                hc = lambda_cond/Dc(L-i)*0.027*Re^0.8*Pr_c^(1/3)
                Tcool(L-i) = Twc(L-i)-q/hc 
                deltaT2 = Tcool(L-i)-Tcool(L+1-i)
            end
            break
        end
    end
    
    if deltaT2 < 0
        Tcool(L-i) = Tcool(L-i+1)
    end
    
    f = 1/(-1.8 * log10(6.9/Re))^2
    ddeltaP = f*0.001/Dc(L-i)*rho_c*uc(L-i).^2/2
    deltaP = deltaP +ddeltaP
    Dp(L-i) = deltaP
    
    if L-i > 1
        Twg(L-1-i) = Twg(L-i)
        Twc(L-1-i) = Twc(L-i)
        Tcool(L-1-i) = Tcool(L-i)
    end
end
max(Tcool)
max(Twg)
max(Dp)