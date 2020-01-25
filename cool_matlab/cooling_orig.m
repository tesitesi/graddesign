% ƒƒPƒbƒgƒGƒ“ƒWƒ“Ä¶—â‹pŒvŽZ

gamma_g = 1.125  %ƒKƒX”ä”M”ä
gamma_c = 1.66      %—â‹pÜ”ä”M”ä
Pr_g = 4 * gamma_g / (9 * gamma_g - 5)    %ƒKƒXƒvƒ‰ƒ“ƒgƒ‹”
Tcin = 112  %—â‹pÜ“üŒû‰·“x [K]
Cpg = 8000  %”RÄƒKƒX’èˆ³”ä”M[J/kg K]
Cpc = 6000   %—â‹pÜ’èˆ³”ä”M [J/kg K]
Rt = 0.05788203 %ƒXƒ[ƒg”¼Œa [m]
Dt = Rt * 2     %ƒXƒ[ƒgŒa [m]
At = Rt^2 * pi  %ƒXƒ[ƒg–ÊÏ
Tc = 3510.96    %ƒ`ƒƒƒ“ƒo‰·“x [K]
g = 9.80655 %d—Í‰Á‘¬“x[m/s^2]
Pc = 4.0 *10^6  %ƒ`ƒFƒ“ƒoˆ³—Í[Pa]
Gc = 6.8182 %—â‹pÜ—¬—Ê[kg/s]
Gg = 30.0   %„iÜ—¬—Ê[kg/s]
c = Pc * At / Gg   %“Á«”r‹C‘¬“x[m/s]
lambda = 390    %”M“`“±—¦[W/mK]
lambda_cond = 0.189  %—â‹pÜ”M’ïR[W/m*K]
mu_c = 0.000123  %—â‹pÜ”S«ŒW”[Pa*s]
rho_c = 426.32
Pr_c = mu_c * Cpc / lambda_cond    %—â‹pÜƒvƒ‰ƒ“ƒgƒ‹”

deltaP = 0

x = cccalcS5(:,1)    %Ž²•ûŒüÀ•W
r = cccalcS5(:,2)   %ƒmƒYƒ‹”¼Œa
A = cccalcS5(:,3)   %ƒmƒYƒ‹–ÊÏ
rho = cccalcS5(:,4) %”RÄƒKƒX–§“x
T = cccalcS5(:,16)   %”RÄƒKƒX‰·“x
p = cccalcS5(:,6)   %”RÄƒKƒXˆ³—Í
M = cccalcS5(:,7)   %”RÄƒKƒXƒ}ƒbƒn”
ug = cccalcS5(:,8)  %”RÄƒKƒX‘¬“x
b = cccalcS5(:,9)   %—â‹pÜ—¬˜H‚‚³
t = cccalcS5(:,10)  %—â‹pÜ—¬˜HŠÔ•
e = cccalcS5(:,11)  %—â‹p•Ç“÷Œú
a = cccalcS5(:,12)  %—â‹pÜ—¬˜H•
Ac = cccalcS5(:,13) %—â‹pÜ—¬˜H–ÊÏ
uc = cccalcS5(:,14)*426/220 %—â‹pÜ—¬‘¬
Dc = cccalcS5(:,15) %—â‹pÜ—¬˜H“™‰¿’¼Œa
lambda = cccalcS5(:,20)
Dp = cccalcS5(:,21)

throat = 745   %ƒXƒ[ƒgˆÊ’u
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