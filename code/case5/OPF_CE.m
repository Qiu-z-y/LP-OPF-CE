function [W1, P1, Q1,Plm1, Qlm1, Vk1, alpha1, W2, P2, Q2,Plm2, Qlm2,Vk2, alpha2, Pdevice, Qdevice, total_cost , loss_cost, device_cost] = myOPF(nsop, nrpfc, fr_sop, to_sop,Cmax_sop, fr_rpfc, to_rpfc, Cmax_rpfc, mpc1,baseMVA1, bus1, gen1, branch1, gencost1, nl1, ns1, ng1, Ybus1, mpc2, baseMVA2, bus2, gen2, branch2, gencost2, nl2, ns2, ng2, Ybus2)
    W1 = sdpvar(2*ns1, 2*ns1,'symmetric');
    W2 = sdpvar(2*ns2, 2*ns2,'symmetric');
    alpha1 = sdpvar(ng1,1);
    alpha2 = sdpvar(ng2,1);
    ndeice = nsop + nrpfc;
    fr = [fr_sop, fr_rpfc];
    to = [to_sop, to_rpfc];
    Cmax = [Cmax_sop, Cmax_rpfc];
    
    
    Pdevice = sdpvar(ndeice,2);
    Qdevice = sdpvar(ndeice,2);
    
   
    
    Pmin_cons1 = [];
    Qmin_cons1 = [];
    Pmax_cons1 = [];
    Qmax_cons1 = [];
    Vmax_cons1 = [];
    Vmin_cons1 = [];
    Pd_cons1 = [];
    Qd_cons1 = [];
    Plm_cons1 = zeros(ns1);
    Slm_cons1 = zeros(ns1);
    
    Pmin_cons2 = [];
    Qmin_cons2 = [];
    Pmax_cons2 = [];
    Qmax_cons2 = [];
    Vmax_cons2 = [];
    Vmin_cons2 = [];
    Pd_cons2 = [];
    Qd_cons2 = [];
    Plm_cons2 = zeros(ns2);
    Slm_cons2 = zeros(ns2);
    
    e1 = eye(ns1);
    e2 = eye(ns2);
    
    P1 = [];
    Q1 = [];
    Plm1 = [];
    Qlm1= [];
    Vk1 = [];
    
    P2 = [];
    Q2 = [];
    Plm2 = [];
    Qlm2= [];
    Vk2 = [];
    
    Closs = 25*baseMVA1;
    
    for i = 1:ns1
        yktmp = e1(:, i)*(e1(:, i))'*Ybus1;
        P1 = [P1;trace(0.5.*[real(yktmp + yktmp.'),  imag(yktmp.' - yktmp); imag(yktmp - yktmp.'), real(yktmp + yktmp.')] * W1)];
        Q1 = [Q1;trace(-0.5.*[imag(yktmp + yktmp.'), real(yktmp - yktmp.'); real(yktmp.' - yktmp), imag(yktmp + yktmp.')] * W1)];
        Vk1 = [Vk1;trace([e1(:, i)*(e1(:,i)'), zeros(ns1); zeros(ns1), e1(:, i)*(e1(:,i)')] * W1)];
    
        Pmin_cons1 = [Pmin_cons1;0];
        Pmax_cons1 = [Pmax_cons1;0];
        Qmin_cons1 = [Qmin_cons1;0];
        Qmax_cons1 = [Qmax_cons1;0];
        Vmin_cons1 = [Vmin_cons1;bus1(i,13)];
        Vmax_cons1 = [Vmax_cons1;bus1(i,12)];
    
        Pd_cons1 = [Pd_cons1;bus1(i,3)];
        Qd_cons1 = [Qd_cons1; bus1(i,4)];
    end
    
    for i = 1:ns2
        yktmp = e2(:, i)*(e2(:, i))'*Ybus2;
        P2 = [P2;trace(0.5.*[real(yktmp + yktmp.'),  imag(yktmp.' - yktmp); imag(yktmp - yktmp.'), real(yktmp + yktmp.')] * W2)];
        Q2 = [Q2;trace(-0.5.*[imag(yktmp + yktmp.'), real(yktmp - yktmp.'); real(yktmp.' - yktmp), imag(yktmp + yktmp.')] * W2)];
        Vk2 = [Vk2;trace([e2(:, i)*(e2(:,i)'), zeros(ns2); zeros(ns2), e2(:, i)*(e2(:,i)')] * W2)];
    
        Pmin_cons2 = [Pmin_cons2;0];
        Pmax_cons2 = [Pmax_cons2;0];
        Qmin_cons2 = [Qmin_cons2;0];
        Qmax_cons2 = [Qmax_cons2;0];
        Vmin_cons2 = [Vmin_cons2;bus2(i,13)];
        Vmax_cons2 = [Vmax_cons2;bus2(i,12)];
    
        Pd_cons2 = [Pd_cons2;bus2(i,3)];
        Qd_cons2 = [Qd_cons2; bus2(i,4)];
    end
    
    
    for i = 1:nv
        Pmax_cons1(busPv(i,1), 1) = Pmax_cons1(busPv(i,1), 1) + Pv(i, 1);
        Pmin_cons1(busPv(i,1), 1) = Pmin_cons1(busPv(i,1), 1) + Pv(i, 2);
        Pmax_cons2(busPv(i,1), 1) = Pmax_cons2(busPv(i,1), 1) + Pv(i, 1);
        Pmin_cons2(busPv(i,1), 1) = Pmin_cons2(busPv(i,1), 1) + Pv(i, 2);
    end
    
    
    for i = 1:ng1
        Qmax_cons1(gen1(i,1),1) = Qmax_cons1(gen1(i,1),1) + gen1(i,4);
        Qmin_cons1(gen1(i,1),1) = Qmin_cons1(gen1(i,1),1) + gen1(i,5);
        Pmax_cons1(gen1(i,1),1) = Pmax_cons1(gen1(i,1),1) + gen1(i,9);
        Pmin_cons1(gen1(i,1),1) = Pmin_cons1(gen1(i,1),1) + gen1(i,10);
    end
    
    for i = 1:ng2
        Qmax_cons2(gen2(i,1),1) = Qmax_cons2(gen2(i,1),1) + gen2(i,4);
        Qmin_cons2(gen2(i,1),1) = Qmin_cons2(gen2(i,1),1) + gen2(i,5);
        Pmax_cons2(gen2(i,1),1) = Pmax_cons2(gen2(i,1),1) + gen2(i,9);
        Pmin_cons2(gen2(i,1),1) = Pmin_cons2(gen2(i,1),1) + gen2(i,10);
    end
    
    for i = 1:nl1
        p = branch1(i,1);
        q = branch1(i,2);
        yltmp = Ybus1(p,q).*(e1(:, p)*(e1(:, p)')-e1(:, p)*(e1(:, q)'));
        Plm1 = [Plm1; trace(0.5.*[real(yltmp + yltmp.'), imag(yltmp.' - yltmp); imag(yltmp - yltmp.'), real(yltmp + yltmp.')] * W1)];
        Qlm1 = [Qlm1; trace(-0.5.*[imag(yltmp + yltmp.'), real(yltmp - yltmp.'); real(yltmp.' - yltmp), imag(yltmp + yltmp.')] * W1)];
        Plm_cons1(p, q) = Plm_cons1(p, q) + branch1(i,6);
        Slm_cons1(p, q) = Slm_cons1(p, q) + branch1(i,6);
    end
    
    for i = 1:nl1
        p = branch2(i,1);
        q = branch2(i,2);
        yltmp = Ybus2(p,q).*(e2(:, p)*(e2(:, p)')-e2(:, p)*(e2(:, q)'));
        Plm2 = [Plm2; trace(0.5.*[real(yltmp + yltmp.'), imag(yltmp.' - yltmp); imag(yltmp - yltmp.'), real(yltmp + yltmp.')] * W2)];
        Qlm2 = [Qlm2; trace(-0.5.*[imag(yltmp + yltmp.'), real(yltmp - yltmp.'); real(yltmp.' - yltmp), imag(yltmp + yltmp.')] * W2)];
        Plm_cons2(p, q) = Plm_cons2(p, q) + branch2(i,6);
        Slm_cons2(p, q) = Slm_cons2(p, q) + branch2(i,6);
    end
    
    
    Plmtmp1 = [];
    Qmltmp1 = [];
    for i = 1:nl1
        p = branch1(i,2);
        q = branch1(i,1);
        yltmp = Ybus1(p,q).*(e1(:, p)*(e1(:, p)')-e1(:, p)*(e1(:, q)'));
        Plmtmp1 = [Plmtmp1; trace(0.5.*[real(yltmp + yltmp.'), imag(yltmp.' - yltmp); imag(yltmp - yltmp.'), real(yltmp + yltmp.')] * W1)];
        Qmltmp1 = [Qmltmp1; trace(-0.5.*[imag(yltmp + yltmp.'), real(yltmp - yltmp.'); real(yltmp.' - yltmp), imag(yltmp + yltmp.')] * W1)];
        Plm_cons1(p, q) = Plm_cons1(p, q) + branch1(i,6);
        Slm_cons1(p, q) = Slm_cons1(p, q) + branch1(i,6);
    end
    Plm1 = [Plm1, Plmtmp1];
    Qlm1 = [Qlm1, Qmltmp1];
    
    Plmtmp2 = [];
    Qmltmp2 = [];
    for i = 1:nl2
        p = branch2(i,2);
        q = branch2(i,1);
        yltmp = Ybus2(p,q).*(e2(:, p)*(e2(:, p)')-e2(:, p)*(e2(:, q)'));
        Plmtmp2 = [Plmtmp2; trace(0.5.*[real(yltmp + yltmp.'), imag(yltmp.' - yltmp); imag(yltmp - yltmp.'), real(yltmp + yltmp.')] * W2)];
        Qmltmp2 = [Qmltmp2; trace(-0.5.*[imag(yltmp + yltmp.'), real(yltmp - yltmp.'); real(yltmp.' - yltmp), imag(yltmp + yltmp.')] * W2)];
        Plm_cons2(p, q) = Plm_cons2(p, q) + branch2(i,6);
        Slm_cons2(p, q) = Slm_cons2(p, q) + branch2(i,6);
    end
    Plm2 = [Plm2, Plmtmp2];
    Qlm2 = [Qlm2, Qmltmp2];
    
    
    Pmin_cons1 = Pmin_cons1 ./ baseMVA1;
    Pmax_cons1 = Pmax_cons1 ./ baseMVA1;
    Qmin_cons1 = Qmin_cons1 ./ baseMVA1;
    Qmax_cons1 = Qmax_cons1 ./ baseMVA1;
    Pd_cons1 = Pd_cons1 ./ baseMVA1;
    Qd_cons1 = Qd_cons1 ./ baseMVA1;
    Plm_cons1 = Plm_cons1 ./ baseMVA1;
    Slm_cons1 = Slm_cons1 ./ baseMVA1;
    
    Pmin_cons2 = Pmin_cons2 ./ baseMVA2;
    Pmax_cons2 = Pmax_cons2 ./ baseMVA2;
    Qmin_cons2 = Qmin_cons2 ./ baseMVA2;
    Qmax_cons2 = Qmax_cons2 ./ baseMVA2;
    Pd_cons2 = Pd_cons2 ./ baseMVA2;
    Qd_cons2 = Qd_cons2 ./ baseMVA2;
    Plm_cons2 = Plm_cons2 ./ baseMVA2;
    Slm_cons2 = Slm_cons2 ./ baseMVA2;
    
    obj = sum(alpha1) + sum(alpha2);
    for i = 1:nv
        p = busPv(i);
        obj = obj + (Pmax_cons1(p) - (Pd_cons1(p) + P1(p,1)))*Cpv;
        obj = obj + (Pmax_cons2(p) - (Pd_cons2(p) + P2(p,1)))*Cpv;
    end
    
    for i = 1:nsop
        obj = obj + 1000*baseMVA1*Cmax(1,i)*1000 / (365*10*24*7);
    end
    for i = 1:nrpfc
        obj = obj + 234*baseMVA1*Cmax(1,i+nsop)*1000 / (365*10*24*7);
    end

    for i = 1:nl1
        obj = obj + Closs * abs(Plm1(i,1)+Plm1(i,2));
    end
    
    for i = 1:nl2
        obj = obj + Closs * abs(Plm2(i,1)+Plm2(i,2));
    end
    
    constraint = [];
    for i = 1:ng1
        ck2 = gencost1(i,5)*baseMVA1^2;
        ck1 = gencost1(i,6)*baseMVA1;
        ck0 = gencost1(i,7);
        alpk = alpha1(i,1);
        p = gen1(i,1);
        F = [ck1*P1(p,1)-alpk+ck0+ck1*Pd_cons1(p,1),  sqrt(ck2)*P1(p,1)+sqrt(ck2)*Pd_cons1(p,1);  sqrt(ck2)*P1(p,1)+sqrt(ck2)*Pd_cons1(p,1), -1];
        constraint = [constraint, F <= 0];
    end
    
    for i = 1:ng2
        ck2 = gencost2(i,5)*baseMVA2^2;
        ck1 = gencost2(i,6)*baseMVA2;
        ck0 = gencost2(i,7);
        alpk = alpha2(i,1);
        p = gen2(i,1);
        F = [ck1*P2(p,1)-alpk+ck0+ck1*Pd_cons2(p,1),  sqrt(ck2)*P2(p,1)+sqrt(ck2)*Pd_cons2(p,1);  sqrt(ck2)*P2(p,1)+sqrt(ck2)*Pd_cons2(p,1), -1];
        constraint = [constraint, F <= 0];
    end
    
    for i = 1:nv
        ck2 = Pvcost(i,1)*baseMVA1^2;
        ck1 = Pvcost(i,2)*baseMVA1;
        ck0 = Pvcost(i,3);
        alpk = alpha1(i+ng1,1);
        p = busPv(i,1);
        F = [ck1*P1(p,1)-alpk+ck0+ck1*Pd_cons1(p,1),  sqrt(ck2)*P1(p,1)+sqrt(ck2)*Pd_cons1(p,1);  sqrt(ck2)*P1(p,1)+sqrt(ck2)*Pd_cons1(p,1), -1];
        constraint = [constraint, F <= 0];
    end
    
    for i = 1:nv
        ck2 = Pvcost(i,1)*baseMVA2^2;
        ck1 = Pvcost(i,2)*baseMVA2;
        ck0 = Pvcost(i,3);
        alpk = alpha2(i+ng2,1);
        p = busPv(i,1);
        F = [ck1*P2(p,1)-alpk+ck0+ck1*Pd_cons2(p,1),  sqrt(ck2)*P2(p,1)+sqrt(ck2)*Pd_cons2(p,1);  sqrt(ck2)*P2(p,1)+sqrt(ck2)*Pd_cons2(p,1), -1];
        constraint = [constraint, F <= 0];
    end
    
    
    for i = 1:nl1
        p = branch1(i,1);
        q = branch1(i,2);
        F = [-Slm_cons1(p,q)*Slm_cons1(p,q), Plm1(i,1), Qlm1(i,1); Plm1(i,1), -1,0; Qlm1(i,1), 0, -1];
        constraint = [constraint, F <= 0];
        constraint = [constraint, Plm1(i,1) <= Plm_cons1(p,q)];
    
        p = branch1(i,2);
        q = branch1(i,1);
        F = [-Slm_cons1(p,q)*Slm_cons1(p,q), Plm1(i,2), Qlm1(i,2); Plm1(i,2), -1,0; Qlm1(i,2), 0, -1];
        constraint = [constraint, F <= 0];
        constraint = [constraint, Plm1(i,2) <= Plm_cons1(p,q)];
    end
    
    for i = 1:nl2
        p = branch2(i,1);
        q = branch2(i,2);
        F = [-Slm_cons2(p,q)*Slm_cons2(p,q), Plm2(i,1), Qlm2(i,1); Plm2(i,1), -1,0; Qlm2(i,1), 0, -1];
        constraint = [constraint, F <= 0];
        constraint = [constraint, Plm2(i,1) <= Plm_cons2(p,q)];
    
        p = branch2(i,2);
        q = branch2(i,1);
        F = [-Slm_cons2(p,q)*Slm_cons2(p,q), Plm2(i,2), Qlm2(i,2); Plm2(i,2), -1,0; Qlm2(i,2), 0, -1];
        constraint = [constraint, F <= 0];
        constraint = [constraint, Plm2(i,2) <= Plm_cons2(p,q)];
    end
    
    j = 1;
    for i = 1:ns1
       if ismember(i,fr)
            constraint = [constraint, Pmin_cons1(i) - Pd_cons1(i) <= P1(i,1) + Pdevice(j,1) <= Pmax_cons1(i) - Pd_cons1(i), Qmin_cons1(i) - Qd_cons1(i) <= Q1(i,1) + Qdevice(j,1) <= Qmax_cons1(i) - Qd_cons1(i)];
            j = j + 1;
       else
            constraint = [constraint, Pmin_cons1(i) - Pd_cons1(i) <= P1(i,1) <= Pmax_cons1(i) - Pd_cons1(i), Qmin_cons1(i) - Qd_cons1(i) <= Q1(i,1) <= Qmax_cons1(i) - Qd_cons1(i)];
       end
       constraint = [constraint, Vmin_cons1(i)*Vmin_cons1(i) <= Vk1(i,1) <= Vmax_cons1(i)*Vmax_cons1(i)];
    end
    constraint = [constraint, W1 >= 0];
    
    
    j = 1;
    for i = 1:ns2
       if ismember(i,to)
            constraint = [constraint, Pmin_cons2(i) - Pd_cons2(i) <= P2(i,1) + Pdevice(j,2) <= Pmax_cons2(i) - Pd_cons2(i), Qmin_cons2(i) - Qd_cons2(i) <= Q2(i,1) + Qdevice(j,2) <= Qmax_cons2(i) - Qd_cons2(i)];
            j = j + 1;
       else
            constraint = [constraint, Pmin_cons2(i) - Pd_cons2(i) <= P2(i,1) <= Pmax_cons2(i) - Pd_cons2(i), Qmin_cons2(i) - Qd_cons2(i) <= Q2(i,1) <= Qmax_cons2(i) - Qd_cons2(i)];
       end
       constraint = [constraint, Vmin_cons2(i)*Vmin_cons2(i) <= Vk2(i,1) <= Vmax_cons2(i)*Vmax_cons2(i)];
    end
    constraint = [constraint, W2 >= 0];
    
    
    for i = 1:ndeice
        if i <= nsop
            constraint = [constraint, Pdevice(i,1) + Pdevice(i,2) == 0];
        else
            constraint = [constraint, Pdevice(i,1) + Pdevice(i,2) == 0, Qdevice(i,1)+Qdevice(i,2) == 0];
        end
        constraint = [constraint, [-Cmax(1,i)*Cmax(1,i), Pdevice(i,1), Qdevice(i,1); Pdevice(i,1), -1, 0; Qdevice(i,1), 0, -1] <= 0];
        constraint = [constraint, [-Cmax(1,i)*Cmax(1,i), Pdevice(i,2), Qdevice(i,2); Pdevice(i,2), -1, 0; Qdevice(i,2), 0, -1] <= 0];
    end
   
    % The constraints and objective function have both been constructed. 
    % At this point, the solver is called and the constraints and objective function are passed into it.
end
