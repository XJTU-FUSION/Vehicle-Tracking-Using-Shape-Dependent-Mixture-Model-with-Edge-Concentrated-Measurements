function [xk_o,Pk_o,betak_o,alphak_o,phik_o,phi_sigmak_o,dirkL_o,C_betak_o,C_alphak_o] = ...
    Update_ad_vel_Angle(xkk_1i,Pkk_1i,betakk_1i,alphakk_1i,C_betakk_1i,C_alphakk_1i,phikk_1i,...
    phi_sigmakk_1i,dirkkL_1i,yk,H,R,L,h1t,h2t,y_velocity,Ns,var_velocity,...
    Ch1t,Ch2t)
% Function: Update_ad_vel_Angle  
% This function performs the measurement update step of the EOT-SD-V
% algorithm using variational Bayesian.

Nmeas = size(yk,2);
scale = 1/3;

G1 = [0,0,1,0;0,0,0,1]';
G2 = [0,1,0,0;-1,0,0,0];

xs = Ns(1,1);
ys = Ns(2,1);
Ms = [-ys;xs];
var_vel = var_velocity; 


Pkk_1 = Pkk_1i; 
xkk_1 = xkk_1i; 
iPkk_1 = Pkk_1^(-1);


phikk_1 = phikk_1i;
phi_sigmakk_1 = phi_sigmakk_1i;


betakk_1 = betakk_1i;
alphakk_1 = alphakk_1i;

C_betakk_1 = C_betakk_1i;
C_alphakk_1 = C_alphakk_1i;

dirkL_1 = dirkkL_1i;



%% Initialization
Pk = Pkk_1;
xk = xkk_1;

phik = phikk_1;
phi_sigmak = phi_sigmakk_1;

betak = betakk_1;
alphak = alphakk_1;

C_betak = C_betakk_1;
C_alphak = C_alphakk_1;

Pzk = cell(1,Nmeas);
zk = yk; 
for jj = 1:Nmeas
    Pzk{1,jj} = R;
end

dirkL = dirkL_1;

rklj_para = zeros(L,Nmeas);
sum_dir = sum(dirkL);
for l = 1:L
    for jj = 1:Nmeas
        rklj_para(l,jj)=dirkL_1(1,l)/sum_dir;
    end
end

iR = R^(-1);
NumOfIterations = 10;
ny = 2;

%% Variational Iterations
for i=1:NumOfIterations 

    GZ2 = cell(1,Nmeas);
    GZ3 = cell(1,Nmeas);
    GZ4 = cell(1,Nmeas);
    EZZ = cell(1,Nmeas);
    for jj =1:Nmeas
        zk_temp = zk(:,jj);
        Pzk_temp = Pzk{1,jj};
        z1 = zk_temp(1);
        z2 = zk_temp(2);
        z1_minus = z1-xs;
        z2_minus = z2-ys;
        mode_sqrt_minus = sqrt((z1_minus)^2+(z2_minus)^2);
        mode_sqrt_3_minus = (mode_sqrt_minus)^3;


        GZ3_temp = [z1_minus/mode_sqrt_minus;z2_minus/mode_sqrt_minus];
        GZ3{1,jj} = GZ3_temp;

        GZ4_temp(1,1) = ((z2_minus)^2)/mode_sqrt_3_minus;
        GZ4_temp(1,2) = (-(z1_minus*z2_minus))/mode_sqrt_3_minus;
        GZ4_temp(2,1) = GZ4_temp(1,2);
        GZ4_temp(2,2) = ((z1_minus)^2)/mode_sqrt_3_minus;
        GZ4{1,jj} = GZ4_temp;

        GZ2_temp = (GZ3_temp-GZ4_temp*zk_temp)';
        GZ2{1,jj} = GZ2_temp;

        EZZ_temp = GZ4_temp*(zk_temp*zk_temp'+Pzk_temp)*GZ4_temp'+...
            GZ4_temp*zk_temp*GZ2_temp+GZ2_temp'*zk_temp'*GZ4_temp'+(GZ2_temp'*GZ2_temp);

        EZZ{1,jj} = EZZ_temp;
    end

    Thetak = phik(1,1);
    SigmaThetak = phi_sigmak(1,1);
    omegak = phik(2,1);
    omega_sigmak = phi_sigmak(2,2);
 
    Rot_matrix = @(tk)[cos(tk) -sin(tk);
        sin(tk)  cos(tk)]; 
    Rot_matrix_diff = @(tk)[-sin(tk) -cos(tk);
        cos(tk)  -sin(tk)];  
 
    [Ck_bar] = cal_truncated_inv_gamma_new(C_alphak,C_betak,1);
    [sqrt_Ck_bar] = cal_truncated_inv_gamma_new(C_alphak,C_betak,0.5);
    [Ck_inv_bar] = cal_truncated_inv_gamma_new(C_alphak,C_betak,-1);
    [sqrt_Ck_inv_bar] = cal_truncated_inv_gamma_new(C_alphak,C_betak,-0.5);


    %% Cht_bar,Cht_inv_bar
    scale_inv = scale^(-1);
    HC1 = [1,0;0,0];HC2 = [0,0;0,1];
    Cht_bar = cell(1,L);
    Cht_inv_bar = cell(1,L);

    for l =1:L
        if mod(l,2)==1 
            Cht_bar{1,l} = scale*diag([1,0])*Ck_bar+scale*diag([0,1]);
            Cht_inv_bar{1,l} = scale_inv*diag([1,0])*Ck_inv_bar+ scale_inv*diag([0,1]);
        elseif mod(l,2)==0 
            Cht_bar{1,l} = scale*diag([0,1])*Ck_bar+scale*diag([1,0]);
            Cht_inv_bar{1,l} = scale_inv*diag([0,1])*Ck_inv_bar+ scale_inv*diag([1,0]);
        end
    end
    %% ln_Cht_bar
    ln_Cht_bar = cell(1,L);

    for l = 1:L
        switch l
            case {1,3}
                ln_Cht_bar{1,l} = cal_ln_Cht(C_alphak(1),C_betak(1),scale);
            case {2,4}
                ln_Cht_bar{1,l} = cal_ln_Cht(C_alphak(2),C_betak(2),scale);
        end
    end

    %% ht_bar
    ht_bar = cell(1,L);
    for l=1:L
        h1t_temp = h1t(:,l);
        h2t_temp = h2t(:,l);
        ht_bar{1,l} = h1t_temp+sqrt_Ck_bar*h2t_temp;
    end

    Cht_inv_ht_bar = cell(1,L);
    scale_temp = scale^(-2);
    for l = 1:L
        Ch1t_temp = Ch1t{1,l};
        Ch2t_temp = Ch2t{1,l};
        h1t_temp = h1t(:,l);
        h2t_temp = h2t(:,l);
        Cht_inv_ht_bar{1,l} = scale_temp*(Ch1t_temp*Ck_inv_bar*h1t_temp+...
            Ch1t_temp*sqrt_Ck_inv_bar*h2t_temp+Ch2t_temp*h1t_temp+Ch2t_temp*sqrt_Ck_bar*h2t_temp);
    end



    %% ht_app
    ht_app = cell(1,L);
    for l = 1:L
        Cht_bar_temp = Cht_bar{1,l};
        Cht_inv_ht_bar_temp = Cht_inv_ht_bar{1,l};
        ht_app{1,l} = Cht_bar_temp*Cht_inv_ht_bar_temp;
    end

    %% h_h_trans_Ch
    h_h_trans_Ch = cell(1,L);
    for l = 1:L
        Ch1t_temp = Ch1t{1,l};
        Ch2t_temp = Ch2t{1,l};
        h1t_temp = h1t(:,l);
        h2t_temp = h2t(:,l);

        part1 = h1t_temp*h1t_temp'*Ch2t_temp+h2t_temp*h2t_temp'*Ch1t_temp;
        part2 = (h1t_temp*h2t_temp'*Ch1t_temp+h2t_temp*h1t_temp'*Ch1t_temp)*sqrt_Ck_inv_bar;
        part3 = (h1t_temp*h2t_temp'*Ch2t_temp+h2t_temp*h1t_temp'*Ch2t_temp)*sqrt_Ck_bar;
        part4 = h1t_temp*h1t_temp'*Ch1t_temp*Ck_inv_bar+h2t_temp*h2t_temp'*Ch2t_temp*Ck_bar;
        part5 = scale_temp*(part1+part2+part3+part4);
        h_h_trans_Ch{1,l} = part5;
    end

    %% 
    % E{x}
    xbar = xk;
    % E{z}
    zjbar = zk;
    % E{theta}
    Thetabar = Thetak;
    % E{X}
    alphak_new = alphak-1;
    Xbar=diag(betak./(alphak_new));

    % sqrt_Xbar
    [axesbar] = cal_Xk_moment(alphak,betak,0.5);
    Xinvbar = diag(alphak./(betak));  

    % sqrt_Xinvbar
    [axesinvbar] = cal_Xk_moment(alphak,betak,-0.5);
    
    % ChtXinvbar 
    ChtXinvbar = cell(1,L);
    for l = 1:L
        Cht_inv_bar_temp = Cht_inv_bar{1,l};
        ChtXinvbar{1,l} = Xinvbar*Cht_inv_bar_temp;
    end
    
    % lnXbar,E{lnX} 
    lnXbar = 0;
    for n = 1:ny
        temp = log(betak(n))-psi(alphak(n));
        lnXbar = lnXbar+temp;
    end

    % lnpibar,E{lnpi}
    lnpibar = zeros(1,L);
    sum_dirkL = sum(dirkL);
    for l = 1:L
        lnpibar(1,l) = psi(dirkL(1,l))-psi(sum_dirkL);
    end

    % Akl,E{L(theta)*ChX^-1*L(theta)^T}
    Akl = cell(1,L);
    for l = 1:L
        Xinvbar_temp = ChtXinvbar{1,l};
        Akl{1,l}=ComputeLXL(Thetabar,SigmaThetak, Xinvbar_temp);
    end

    % ukl
    ukl = cell(1,L);
    Rot_temp = ComputeTbar(Thetabar,SigmaThetak);
    for l = 1:L
        temp1 = (Akl{1,l})^(-1);
        temp4 = Cht_inv_ht_bar{1,l};
        ukl{1,l} = temp1*Rot_temp*axesinvbar*temp4;
    end

    % jiankl
    jiankl = cell(1,L);
    for l = 1:L
        temp0 = zeros(2,2);
        for jj = 1:Nmeas
            temp1 = rklj_para(l,jj)*Akl{1,l};
            temp0 = temp0+temp1;
        end
        jiankl{1,l} = temp0;
    end

    % jiankl_sum
    jiankl_sum = [];
    for l = 1:L
        jiankl_temp = jiankl{1,l};
        jiankl_sum = blkdiag(jiankl_sum,jiankl_temp);
    end

    % ukl_bend
    ukl_bend = cell(1,L);
    empty_L = [];
    for l = 1:L
        part1 = zeros(2,1);
        part2 = 0;
        for jj = 1:Nmeas
            part1 = part1 + rklj_para(l,jj)*zk(:,jj);
            part2 = part2 + rklj_para(l,jj);
        end
        if part2==0
            empty_L = [empty_L,l];
            ukl_bend{1,l} = [];
        else
            ukl_bend{1,l} = (part1/part2)-ukl{1,l};
        end
    end

    % ukl_bend_sum
    ukl_bend_sum = [];
    for l = 1:L
        if (~isempty(ukl_bend{1,l}))
            ukl_bend_temp = ukl_bend{1,l};
            ukl_bend_sum =[ukl_bend_sum;ukl_bend_temp];
        end
    end
    
    % HL_sum_ukl
    HL_sum_ukl = [];
    for l = 1:L
        if (~isempty(ukl_bend{1,l}))
            HL_temp = H;
            HL_sum_ukl =[HL_sum_ukl;HL_temp];
        end
    end

    % HL_sum
    HL_sum = [];
    for l = 1:L
        HL_temp = H;
        HL_sum =[HL_sum;HL_temp];
    end

    % jiankl_sum_ukl
    jiankl_sum_ukl = [];
    for l = 1:L
        index_l = find(empty_L==l);
        if isempty(index_l)
            jiankl_temp = jiankl{1,l};
            jiankl_sum_ukl = blkdiag(jiankl_sum_ukl,jiankl_temp);
        end
    end

    %% Update q_rk
    rklj_para_out = zeros(L,Nmeas);
    rklj_para_out_exp = zeros(L,Nmeas);
    rklj_para_out_noexp = zeros(L,Nmeas);
    for jj = 1:Nmeas
        for l = 1:L
            ln_Cht_bar_temp = ln_Cht_bar{1,l};
            part1 = lnpibar(1,l);
            part2 = -0.5*ln_Cht_bar_temp;

            part4 = zjbar(:,jj)-H*xbar-ukl{1,l};
            part5 = Akl{1,l};
            part6 = part4'*part5*part4; 
            part7_new = trace(h_h_trans_Ch{1,l});

            part8 = Pzk{1,jj} + H*Pk*H';
            part9 = trace(part5*part8);
            temp1 = ukl{1,l};
            part10_later = -temp1'*part5*temp1;
            part10 = -0.5*(part6+part10_later+part7_new+part9);
            part11 = part1+part2+part10;
            rklj_para_out_noexp(l,jj) = part11;
            rklj_para_out_exp(l,jj) = exp(part11);
        end
    end


    rklj_para_out_next = rklj_para_out_exp;

    for jj = 1:Nmeas
        sum_jj = sum(rklj_para_out_next(:,jj));
        for l = 1:L
            rklj_para_out(l,jj) = rklj_para_out_next(l,jj)./sum_jj;
        end
    end
    rklj_para_out_4xiaoshu = roundn(rklj_para_out,-4);
    rklj_para_out = rklj_para_out_4xiaoshu;
    rklj_para = rklj_para_out;
    
    % jiankl
    jiankl = cell(1,L);
    for l = 1:L
        temp0 = zeros(2,2);
        for jj = 1:Nmeas
            temp1 = rklj_para(l,jj)*Akl{1,l};
            temp0 = temp0+temp1;
        end
        jiankl{1,l} = temp0;
    end

    % jiankl_sum
    jiankl_sum = [];
    for l = 1:L
        jiankl_temp = jiankl{1,l};
        jiankl_sum = blkdiag(jiankl_sum,jiankl_temp);
    end

    % ukl_bend
    ukl_bend = cell(1,L);
    empty_L = [];
    for l = 1:L
        part1 = zeros(2,1);
        part2 = 0;
        for jj = 1:Nmeas
            part1 = part1 + rklj_para(l,jj)*zk(:,jj);
            part2 = part2 + rklj_para(l,jj);
        end
        if part2==0
            empty_L = [empty_L,l];
            ukl_bend{1,l} = [];
        else
            ukl_bend{1,l} = (part1/part2)-ukl{1,l};
        end
    end

    % ukl_bend_sum
    ukl_bend_sum = [];
    for l = 1:L
        if (~isempty(ukl_bend{1,l}))
            ukl_bend_temp = ukl_bend{1,l};
            ukl_bend_sum =[ukl_bend_sum;ukl_bend_temp];
        end
    end
    % HL_sum_ukl
    HL_sum_ukl = [];
    for l = 1:L
        if (~isempty(ukl_bend{1,l}))
            HL_temp = H;
            HL_sum_ukl =[HL_sum_ukl;HL_temp];
        end
    end

    % HL_sum
    HL_sum = [];
    for l = 1:L
        HL_temp = H;
        HL_sum =[HL_sum;HL_temp];
    end

    % jiankl_sum_ukl
    jiankl_sum_ukl = [];
    for l = 1:L
        index_l = find(empty_L==l);
        if isempty(index_l)
            jiankl_temp = jiankl{1,l};
            jiankl_sum_ukl = blkdiag(jiankl_sum_ukl,jiankl_temp);
        end
    end

    % Gx
    Gx = xk*xk'+Pk;

    % B_j
    B_j = cell(1,Nmeas);
    for jj = 1:Nmeas
        tempzx = zjbar(:,jj)-H*xbar;

        M_inv = tempzx*tempzx'+ Pzk{1,jj} + H*Pk*H';
        B_j{1,jj} = ComputeLXL(-Thetak,SigmaThetak,M_inv);
    end

    Tbar=ComputeTbar(Thetabar,SigmaThetak);
    Ttransbar=Tbar';


    %%  Update q_Ck
    % C_alpha
    temp_rjl_ck1 = 0;
    temp_rjl_ck2 = 0;
    C_alphakk_1_ck1 = C_alphakk_1(1,1);
    C_alphakk_1_ck2 = C_alphakk_1(2,1);
    for l =1:4 
        if mod(l,2)==1 
            for jj = 1:Nmeas
                rklj_para_temp = rklj_para(l,jj);
                temp_rjl_ck1 = temp_rjl_ck1 + rklj_para_temp;
            end
        elseif mod(l,2)==0 
            for jj = 1:Nmeas
                rklj_para_temp = rklj_para(l,jj);
                temp_rjl_ck2 = temp_rjl_ck2 + rklj_para_temp;
            end
        end

    end

    C_alphak_out_ck1 = C_alphakk_1_ck1 + 0.5*temp_rjl_ck1;
    C_alphak_out_ck2 = C_alphakk_1_ck2 + 0.5*temp_rjl_ck2;

    C_alphak_out = [C_alphak_out_ck1;C_alphak_out_ck2];

    % C_beta
    C_gama_beta_ck1 = 0;
    C_gama_beta_ck2 = 0;
    
    for l = 1:4
        for jj  = 1:Nmeas
            rklj_temp= rklj_para(l,jj);%

            part0 = rklj_temp;

            part1 = Xinvbar*B_j{1,jj};

            z_hx_temp = zjbar(:,jj)-H*xbar;
           
            ht_bar_temp = ht_app{1,l};

            part2 = -axesinvbar*Ttransbar*z_hx_temp*(ht_bar_temp');

            part3 = -axesinvbar*ht_bar_temp*(z_hx_temp')*Tbar;

            part4 = axesinvbar*ht_bar_temp*(ht_bar_temp')*axesbar;
            
            if mod(l,2)==1 
                part_sum_ck1 = (1/scale)*part0*(part1+part2+part3+part4)*HC1;
                C_gama_beta_ck1 = C_gama_beta_ck1+trace(part_sum_ck1);
                
            elseif mod(l,2)==0 
                part_sum_ck2 = (1/scale)*part0*(part1+part2+part3+part4)*HC2;
                C_gama_beta_ck2 = C_gama_beta_ck2+trace(part_sum_ck2);
                
            end
        end
    end

    C_betak_out(1,1) = C_betakk_1(1,1) + 0.5*C_gama_beta_ck1;
    C_betak_out(2,1) = C_betakk_1(2,1) + 0.5*C_gama_beta_ck2;

    

   

    %% Update xk
    [xk_vel_cov,xk_vel_mean] =...
        VB_vel_xk(EZZ,G1,G2,GZ2,GZ4,zk,y_velocity,var_vel,Ms,Nmeas,omegak,omega_sigmak);

    Pk_out = (iPkk_1 + HL_sum'*jiankl_sum*HL_sum+xk_vel_cov)^(-1);
    xk_out = Pk_out*(iPkk_1*xkk_1+ HL_sum_ukl'*jiankl_sum_ukl*ukl_bend_sum+xk_vel_mean);

    %% Update zk
    zk_out = zeros(2,Nmeas);
    Pzk_out = cell(1,Nmeas);
    
    part2_E1 = G1'*Gx*G2'+G2*Gx*G1+G1'*xk*Ms'+Ms*xk'*G1;
    part3_E1 = omegak^2+omega_sigmak;
    part4_E1 = G2*Gx*G2'+G2*xk*Ms'+Ms*xk'*G2'+Ms*Ms';
    E1 = G1'*Gx*G1 + omegak*part2_E1+part3_E1*part4_E1;
    E2 = (G1'+omegak*G2)*xk+omegak*Ms;

    for jj = 1:Nmeas
        temp1 = zeros(2,2);
        temp2 = zeros(2,1);
        for l = 1:L
            temp0 = rklj_para(l,jj)*Akl{1,l};
            temp1 = temp1+temp0;
            temp3 = H*xk+ukl{1,l};
            temp2 = temp2+temp0*temp3;
        end

        % covariance
        GZ4_temp = GZ4{1,jj};
        temp_cov = (1/(var_vel))*GZ4_temp'*E1*GZ4_temp;

        % mean
        GZ4_temp = GZ4{1,jj};
        GZ2_temp = GZ2{1,jj};
        mean_temp1 = (-1/(2*var_vel))*(GZ2_temp*E1*GZ4_temp)';
        mean_temp2 = (-1/(2*var_vel))*(GZ4_temp'*E1*GZ2_temp'); %
        mean_temp3 = (y_velocity(jj)/(var_vel))*GZ4_temp'*E2;
        temp_mean = mean_temp1+mean_temp2+mean_temp3;


        Pzk_out{1,jj} = (iR+temp1+temp_cov)^(-1);
        ykj = yk(:,jj);
        zk_out(:,jj) = Pzk_out{1,jj}*(iR*ykj+temp2+temp_mean);
    end


    %% Update pik
    dirkL_out = zeros(1,L);
    for l =1:L
        temp1 = 0;
        for jj = 1:Nmeas
            rklj_para_temp = rklj_para(l,jj);
            temp1 = temp1 + rklj_para_temp;
        end
        dirkL_out(1,l) = dirkL_1(1,l) + temp1;
    end

    %% Update Xk 
    temp_rjl = 0;
    for l =1:L
        for jj = 1:Nmeas
            rklj_para_temp = rklj_para(l,jj);
            temp_rjl = temp_rjl + rklj_para_temp;
        end
    end
    alphak_out = alphakk_1 + 0.5*temp_rjl;
    axes_constant = Xbar*axesinvbar;

    gama_beta = zeros(2,2);
    for l = 1:L
        for jj  = 1:Nmeas
            temp2= rklj_para(l,jj);
            temp4 = zjbar(:,jj)-H*xbar;
            part1 = Cht_inv_bar{1,l}*B_j{1,jj};
            part2 = axes_constant*(h_h_trans_Ch{1,l})*axes_constant';
            part3 = -Ttransbar*temp4*Cht_inv_ht_bar{1,l}'*axes_constant';
            part4 = -axes_constant*Cht_inv_ht_bar{1,l}*temp4'*Tbar;
            gama_beta = gama_beta+temp2*(part1+part2+part3+part4);

        end
    end

    betak_out(1,1) = betakk_1(1,1) + 0.5*gama_beta(1,1);
    betak_out(2,1) = betakk_1(2,1) + 0.5*gama_beta(2,2);

    %% Update phik
    % thetak
    Wlj = cell(1,Nmeas);
    for jj = 1:Nmeas
        inno_temp = zk(:,jj)-H*xk;
        temp_hph = H*Pk*H';
        Wlj{1,jj} = inno_temp*inno_temp'+Pzk{1,jj}+temp_hph;
    end
    innotemp = 0;
    for jj = 1:Nmeas
        for l = 1:L
            %rBk
            part01 = rklj_para(l,jj);
            part02 = ChtXinvbar{1,l}; 
            temp_rBk = part01*part02;

            % Eab1

            part11 = (Rot_matrix_diff(Thetabar))';
            part12 = Wlj{1,jj};
            part13 = Rot_matrix_diff(Thetabar);
            part14 = Thetabar;
            temp_sum1 = part11*part12*part13*part14;

            % Eab2

            part21 = (ChtXinvbar{1,l})^(-1);
            part22 = axesinvbar;
            part23_24_new = Cht_inv_ht_bar{1,l};

            part25 = (zk(:,jj)-H*xk)'*Rot_matrix_diff(Thetabar);
            temp_sum2 = part21*part22*part23_24_new*part25;

            % Eab3
            part31 = (Rot_matrix(Thetabar))';
            temp_sum3 = part31*part12*part13;
            
            sum_ab = temp_sum1 + temp_sum2 - temp_sum3;
            m = trace(temp_rBk*sum_ab);
            innotemp = innotemp + m;
            
        end
    end

    Sigmatemp = 0;
    for jj = 1:Nmeas
        for l = 1:L
            %rBk
            part01 = rklj_para(l,jj);
            part02 = ChtXinvbar{1,l};
            temp_rBk = part01*part02;

            % Ebb1
            part11 = (Rot_matrix_diff(Thetabar))';
            part12 = Wlj{1,jj};
            part13 = Rot_matrix_diff(Thetabar);
            sum_bb = part11*part12*part13;
            Sigmatemp = Sigmatemp + trace(temp_rBk*sum_bb);
        end
    end
   
    % omegak
    omega_sigmatemp = 0;
    omega_meantemp = 0;
    for jj = 1:Nmeas
        zk_temp = zk(:,jj);
        % covariance
        sigma_temp1 = G2*Gx*G2'+G2*xk*Ms'+Ms*xk'*G2'+Ms*Ms';
        EZZ_temp = EZZ{1,jj};
        sigma_temp2 = (1/(var_vel))*trace(sigma_temp1*EZZ_temp);
        omega_sigmatemp = omega_sigmatemp+sigma_temp2;

        % mean
        GZ4_temp = GZ4{1,jj};
        GZ2_temp = GZ2{1,jj};

        mean_temp1 = (y_velocity(1,jj)/var_vel)*(zk_temp'*GZ4_temp'+GZ2_temp)*(G2*xk+Ms);
        mean_temp2 = (-1/var_vel)*trace(G2'*EZZ_temp*G1'*Gx);  
        mean_temp3 = (-1/var_vel)*(Ms'*EZZ_temp*G1'*xk);
        mean_temp4 = mean_temp1+mean_temp2+mean_temp3;

        omega_meantemp = omega_meantemp+mean_temp4;
    end

    cov_phi = [Sigmatemp,0;0,omega_sigmatemp];
    mean_phi = [innotemp;omega_meantemp];
    phi_sigmak_out = (phi_sigmakk_1^(-1)+ cov_phi)^(-1);
    phik_out= phi_sigmak_out*((phi_sigmakk_1^(-1))*phikk_1+mean_phi);
    

    %% Iteration
    zk = zk_out;
    Pzk = Pzk_out;
    
    betak = betak_out;
    alphak = alphak_out;

    C_betak = C_betak_out;
    C_alphak = C_alphak_out;

    dirkL = dirkL_out;

    Pk = Pk_out; 
    xk = xk_out; 
 
    rklj_para = rklj_para_out;

    phi_sigmak = phi_sigmak_out;
    phik= phik_out;
end

%% Output the states
xk_o = xk;
Pk_o = Pk;

betak_o = betak;
alphak_o = alphak;

C_betak_o = C_betak;
C_alphak_o = C_alphak;

dirkL_o = dirkL;
phi_sigmak_o = phi_sigmak;
phik_o = phik_out;
end

function LXLbar=ComputeLXL(Thetabar,SigmaThetak,Xinvbar) 
% Compute E{L(theta)*Xinvbar*L(theta)^T}
m11 = Xinvbar(1,1); m12 = Xinvbar(1,2); m21 = Xinvbar(2,1);  m22 = Xinvbar(2,2);
k11 = (1+cos(2*Thetabar)*exp(-2*SigmaThetak));
k12 = (1-cos(2*Thetabar)*exp(-2*SigmaThetak));
k13 = (sin(2*Thetabar)*exp(-2*SigmaThetak));

LXLbar(1,1) = m11*k11 + m22*k12 - (m21+m12)*k13;
LXLbar(1,2) = m12*k11 - m21*k12 + (m11-m22)*k13;
LXLbar(2,1) = m21*k11 - m12*k12 + (m11-m22)*k13;
LXLbar(2,2) = m22*k11 + m11*k12 + (m21+m12)*k13;

LXLbar=0.5.*LXLbar;

end


function Tbar=ComputeTbar(Thetabar,SigmaThetak)
% Compute the expectation of the rotation matrix
k11 = cos(Thetabar)*exp(-0.5*SigmaThetak);
k12 = -sin(Thetabar)*exp(-0.5*SigmaThetak);
k21 = -k12;
k22 = k11;
Tbar(1,1) = k11;
Tbar(1,2) = k12;
Tbar(2,1) = k21;
Tbar(2,2) = k22;

end