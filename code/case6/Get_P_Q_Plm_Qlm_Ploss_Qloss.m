function [P, Q, P_lm, Q_lm, P_loss, Q_loss] = power_flow(baseMVA, bus, gen, branch, Ybus, V_pu, theta_rad)
    n = length(V_pu);
    V_complex = V_pu .* exp(1j * theta_rad);
    S_inj_pu = V_complex .* conj(Ybus * V_complex);
    
    P = real(S_inj_pu) * baseMVA;
    Q = imag(S_inj_pu) * baseMVA;
    
    nl = size(branch, 1);
    S_fr = zeros(nl, 1);
    S_to = zeros(nl, 1);
    loss = zeros(nl, 1);
    
    for k = 1:nl

        if branch(k, 11) == 0
            S_fr(k) = 0;
            S_to(k) = 0;
            loss(k) = 0;
            continue;
        end
        
        f = branch(k, 1);
        t = branch(k, 2);
        R = branch(k, 3);
        X = branch(k, 4);
        B = branch(k, 5);
        tap = branch(k, 9);
        shift_deg = branch(k, 10);
        
        if tap == 0
            tap = 1;
        end
        shift = shift_deg * pi/180;
        t_ratio = tap * exp(1j * shift);
        

        y_series = 1 / (R + 1j * X);
        y_shunt = 1j * B / 2;
        If = y_series * (V_complex(f) - V_complex(t)/t_ratio)/conj(t_ratio) + V_complex(f) * y_shunt;
        S_fr_pu = V_complex(f) * conj(If);
        It = y_series * (V_complex(t) - V_complex(f)*conj(t_ratio)) + V_complex(t) * y_shunt;
        S_to_pu = V_complex(t) * conj(It);
        S_fr(k) = S_fr_pu * baseMVA;
        S_to(k) = S_to_pu * baseMVA;
        
        P_lm = real(S_fr);
        Q_lm = imag(S_fr);
        P_loss = sum(real(S_fr) + real(S_to));
        Q_loss = sum(imag(S_fr) + imag(S_to));
    end
end
