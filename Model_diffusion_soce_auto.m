clear all

%geometrical parameters and diffusion constant
R_max=0.5; %radius
L_max=1.25; %length half sarcomere
SR_A = 0.54; %mmicrom square
n=10; %sections in radial direction
m=20; %sections in longitudinal direction
D=3*1e2; %microm^2s-1



% Diffusional rate constants and tim e
slow_fact=1;
K_D_R=slow_fact*D*(n/R_max)^2;
K_D_L=slow_fact*D*(m/L_max)^2;


%%%%%%%%%%%%%%VOLUMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_max=pi*R_max^2*L_max; %Vol micrometer^3
V_TC=0.035*V_max; %Vol TC
V_SR=V_max*0.055/(m-1); %Vol long. SR
V_EF=((n-1-1+0.5)^2-(max(0,n-1-2+0.5))^2)/n^2*V_max/m; %Vol each comparterment in the  extrafilm space
%%%%%Each i-th ring has a volume of Vn-1=2*(i-1)/n^2*pi*Rmax^2; if
%%%%%n=31 V_EF=6% Vmax
V_tub=0.01*V_max; %Vol tubule

pos_SR=0.5; %triad at 500 nm from the Z-line (Brown et al., 1998)(Gomez et al. (2006)
comp_SR=round(pos_SR/(L_max/m));
m_SR=m-comp_SR; %position of the RyR

l_actin=1; %Actin length 1 µm
ncomp_Tr=round(l_actin/(L_max/m)); %Number of m compartments which include troponin
m_Tr=m-ncomp_Tr;%number of m sub-comp without Troponin close to the M-line

for ind_n=1:(n-1)
    for ind_m=1:m
        volumes(ind_m+m*(ind_n-1))=((ind_n-1+0.5)^2-(max(0,ind_n-2+0.5))^2)/n^2*V_max/m;
    end
end
V_cito=sum(volumes(1:m*(n-1)));
volumes((m*(n-1))+1:m*(n-1)+m_SR-1)=V_SR;
volumes(m*(n-1)+m_SR)=V_TC;
volumes((m*(n-1)+m_SR+1):(m*n))=V_SR;






p_soce=1;
mcu_syl=1;


for anls=1:2 %with or without PVA
    hertz=[60];% 20 5 1]; %Hz of stimulation
    held=[2];% 5 10 15]; %seconds of stimulation
    if anls == 1
        anl_p=0;
        pvb_p=1;
    else
        pvb_p=0;
        anl_p=0;
    end
    
    if anls == 1
        V_mito=0.042*V_max;
    elseif anls == 2
        V_mito=0.056*V_max; %Feliciano data
    end
    
    Ca_SR_res=0.5*1e3; %micromolar resting concentration in SR
    Ca_myo_res=0.1; %micromolar resting concentration in myoplasm citosol
    
    
    index=1:m*n;
    M_index=reshape(index,m,n);
    M_PVMG_internal_index=reshape(M_index(2:m-1,2:n-2),1,[]);
    M_Tr_index=reshape(M_index(m_Tr+1:m,1:n-2),1,[]);
    
    V_Tr=0;%Volum occupied by Troponin
    for ind_n=1:(n-2)
        for ind_m=1:m
            V_Tr=V_Tr+(ind_m>m_Tr)*((ind_n-1+0.5)^2-(max(0,ind_n-2+0.5))^2)/n^2*V_max/m;
        end
    end
    
    Voltofunc=struct('v', {V_max, V_Tr, V_TC, V_SR, V_EF, V_mito, V_tub});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%Parameters for Citosolic Buffers%%%%%%%%%%%%%%%%%%%%%%%%
    
    TtubB_tot = 400*V_max/V_tub; %Ttubule buffering total capacity Owen JPhys 1997
    n_BTtub = 1; %Ttubule buffering cooperativity
    Fura_tot = 2*V_max/V_cito; %citosolic µM concentration of Fura-2
    
    CaP_fact=pvb_p*1; %Parvalbumin for Calcium
    CaMg_fact=CaP_fact*1; %Parvalbumin for Mg
    CaFura_fact = 1;
    CaT_fact=1;
    
    Q10tr=2^((25-20)/10);
    Q10pv=2^((25-20)/10);
    n_trp1=1; %cooperativity troponin
    n_trp2=1; %cooperativity troponin
    trp1_ratio=2040/500;%8.72;%µM ratio off/on troponin1 from Baylor 2007 JGP
    trp2_ratio=22/115;%0.194;%µM ratio off/on troponin1 from Baylor 2007 JGP
    
    P_tot=1500;%1000; %Total Parv microM BH2007
    Tr_tot=120*V_max/V_Tr;% %Total Troponin microM
    K_D_R_Mg = K_D_R; %Mg radial diffusion
    K_D_L_Mg = K_D_L; %Mg longitudinal diffusion
    ktr1_on = 500*Q10tr;%CaT_fact*177;% microM-1 s-1from Baylor 2007 JGP
    ktr1_off= ktr1_on*trp1_ratio;%CaT_fact*115;%120; %Ca_T off s-1
    ktr2_on = 115*Q10tr;%Barclay 2021, 1.41 fattore Q10 da 20 a 25 gradi
    ktr2_off= ktr2_on*trp2_ratio;%CaT_fact*115;%120; %Ca_T off s-1
    kp_on = CaP_fact*55*Q10pv;%Barclay 2020 0.417*1e2;%*4.5*1e2;%*0.417*1e2;%2.5e2; %Ca_P on microM-1 s-1
    kp_off=CaP_fact*1.2*Q10pv;% Barclay 2020 0.5;%1; %Ca_P off s-1
    kf_on = CaFura_fact*146; % for FURA-2 microM-1 s-1 from Schneider 1988
    kf_off = CaFura_fact*21.9; %for FURA-2 s-1 from Schneider 1988
    km_on=CaMg_fact*0.04*Q10pv;%1.2*1e-1;%6.6e-2; % Mg_P on microM-1 s-1
    km_off=CaMg_fact*4.*Q10pv;%3.4;%3.;%6; % Mg_p off s-1
    k_soce_off = 0.1;
    CaTtub_basal = 1e3; %µM local concentration free calcium in Ttubule basal
    CaBTtub_basal=CaTtub_basal*10; %calcium bound to Ttubule B basal
    kTtb_on =  1.e1; % Ttubule buffering on microM s-1
    kTtb_off = kTtb_on*(CaTtub_basal^n_BTtub*(TtubB_tot-CaBTtub_basal)/CaBTtub_basal);
    
    Bcitotofunc=struct('v', {CaP_fact, CaMg_fact,...
        CaT_fact, P_tot, Fura_tot, kf_on, kf_off, M_PVMG_internal_index, Tr_tot, M_Tr_index, K_D_R_Mg, K_D_L_Mg, ktr1_on, ktr1_off, ktr2_on, ktr2_off, n_trp1, n_trp2, kp_on, kp_off,...
        km_on, km_off, TtubB_tot, n_BTtub, kTtb_on, kTtb_off, k_soce_off});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%Mito buffering definition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%data from figure 2B Bazil et al 2013 %%%%%%%%
    CaMito_basal=0.16; %µM
    MB1_fact=10; %speed factor for first
    MB2_fact=10; %speed factor for second
    n_mb2=6;  %cooperativity second
    MB1_tot=1*35.0*1e3; %Total first
    MB2_tot=1*90.*1e3; %Total second
    Kd_1=2.0; %µM
    Kd_2=1.7; %µM
    
    kmb1_on = MB1_fact*1e1; %  micorM-1 s-1
    kmb1_off = Kd_1*kmb1_on; %s-1
    kmb2_on = MB2_fact*1e1; % micorM-1 s-1
    kmb2_off = Kd_2^n_mb2*kmb2_on; % s-1
    
    MitoBtofunc = struct('v', {MB1_tot, MB2_tot, kmb1_on, kmb1_off, kmb2_on, kmb2_off, n_mb2});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%CSQ definition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CaS2_fact=0.1; %speed factor for second SR buffer
    CaS_fact=1; %speed factor for calsequestrin
    n_csq=3;
    S_tot=82.8*1e3; %Total CalSeq microM =~1.2mM*V_max/V_TC
    ratio_2=1;
    %note: C&A84 S_tot=31mM, BUT  Murphy e Lamb JP 2009 CalS=36
    %umoli/liter of fiber
    %wt 70-80 Ca sites/molecule (Mol Biosyst. 2013 Jul 4; 9(7): 1949?1957.)
    %compared to V_TC means a 83000 for 80 and 72000 for 70
    S2_tot=0.1*82.8*1e3; %Total CalSeq microM =~1.2mM*V_max/V_TC
    ks_on = CaS_fact*3e-3; %Ca_S on  micorM-1 s-1
    ks_off = CaS_fact*3e3; %Ca_Soff s-1
    ks2_on = CaS2_fact*5*3e-3; %Ca_S2 on  micorM-1 s-1
    ks2_off = ks2_on/ratio_2; %Ca_S2off s-1
    
    CSQtofunc = struct('v', {S_tot, S2_tot, ks_on, ...
        ks_off, ks2_on, ks2_off, n_csq, CaS2_fact, CaS_fact});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%Definition of the indexing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mg_b=m*n; %Mg variables starts from m*n +1 till m*n+m*(n-1)
    CaPV_b=Mg_b+m*(n-1); %PV-Ca starts from m*n +1 till m*n+m*(n-1) +Mg_b
    MgPV_b=CaPV_b+m*(n-1); %PV-Mg starts from m*n+m*(n-1)+1 till m*n+2m*(n-1)+Mg_b
    Fura_b=MgPV_b+m*(n-1); %Fura-2 starts from the end of MgPV
    Tr1_b=Fura_b+m*(n-1); %First binding site on Tr starts from m*n+2m*(n-1)+1 till m*n+2m*(n-1) +m*(n-2)+Mg_b
    Tr2_b=Tr1_b+m*(n-2); %Second binding site on Tr starts from m*n+2m*(n-1)+1 till m*n+2m*(n-1) +m*(n-2)+Mg_b
    CaS_b=Tr2_b+m*(n-2); %CaS starts from m*n+2m*(n-1) +m*(n-2)+1+Mg_b
    CaS=CaS_b+1; %Calsequestrin index
    CaS2=CaS_b+2; %Second RS buffer index
    Mito_free = CaS_b+3; %free Ca in Mitochondria Index
    Mito_bound1 = CaS_b+4; %Index bound Ca in inner mitochondrial space
    Mito_bound2 = CaS_b+5; %Index bound Ca in inner mitochondrial space
    Ttubule_free = CaS_b+6; %Index Ttuble calcium free
    Ttubule_bound = CaS_b+7; %Index Ttuble calcium free
    all_var=Ttubule_bound; %all variables index
    
    indextofunc=struct('v', {Mg_b, CaPV_b, MgPV_b, Fura_b, Tr1_b, Tr2_b, CaS_b, CaS, CaS2,...
        Mito_free, Mito_bound1, Mito_bound2, Ttubule_free, Ttubule_bound, all_var});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for freq=1:length(hertz)
        if anls == 1
            filename=sprintf('A_WT_Hz%u_S%u_MCU',hertz(freq),mcu_syl);
        else
            filename=sprintf('A_CSQKO_Hz%u_S%u_MCU',hertz(freq),mcu_syl);
        end
        clear x_s t_s x0
        tic
        
        % Diffusional rate constants and time
        t_start=0;
        t_finish=held(freq)+4; %total time of the analysis
        t_step=0.1; %length of each substeps
        t_ini=0; %initial time
        t_end=t_step; %end time for each substep
        t_sampling=1000; %number of sampling points per second
        
        
        
        ind_cycle=-1; %index of the next while cycle
        avg=0;
        while t_end<=t_finish+t_step/2
            ind_cycle=ind_cycle+1;
            tspan=[t_ini t_end]; %time span
            % tspan = linspace(t_ini, t_end, t_points);
            
            if t_ini==0
                
                x0=zeros(all_var,1); %initial concentration matrix
                
                x0(1:m*(n-1))=Ca_myo_res;
                x0(m*(n-1)+1:m*n)=Ca_SR_res;
                
                Mg_basal=1.e3; % basal mg concentration µM
                x0(Mg_b+1:CaPV_b) = Mg_basal; %Total Magnesium microM
                
                if pvb_p
                    x0(CaPV_b+1:MgPV_b) = Ca_myo_res*kp_on/kp_off*P_tot/(1+Ca_myo_res*kp_on/kp_off+Mg_basal*km_on/km_off); %Initial [CaP] concentration microM in EF
                    x0(MgPV_b+1:Tr1_b) =Mg_basal*km_on/km_off*P_tot/(1+Mg_basal*km_on/km_off+Ca_myo_res*kp_on/kp_off);%Initial [MgP] concentration microM in Myo
                else
                    x0(CaPV_b+1:MgPV_b) = Ca_myo_res*4.5e2*1000/(1+Ca_myo_res*4.5e2+Mg_basal*1.2e-1/3.4); %Initial [CaP] concentration microM in EF
                    x0(MgPV_b+1:Tr1_b) =Mg_basal*1.2e-1/3.4*1000/(1+Mg_basal*1.2e-1/3.4+Ca_myo_res*4.5e2);%Initial [MgP] concentration microM in Myo
                end
                
                x0(Fura_b+1:Tr1_b) = kf_on*Ca_myo_res*Fura_tot/(kf_off+kf_on*Ca_myo_res);
                
                Tr1_eq=Tr_tot/(trp1_ratio/Ca_myo_res^n_trp1+1+Ca_myo_res^n_trp1/trp2_ratio);
                x0(Tr1_b+1:Tr2_b) = Tr1_eq; %Initial [CaT] concentration microM
                
                x0(Tr2_b+1:CaS_b) = Ca_myo_res^n_trp2*Tr1_eq/trp2_ratio; %Initial [CaT] concentration microM
                
                x0(CaS) = ks_on*Ca_SR_res^n_csq*S_tot/(ks_off+ks_on*Ca_SR_res^n_csq); %Calsequestrin
                
                x0(CaS2) = ks2_on*Ca_SR_res*S2_tot/(ks2_off+ks2_on*Ca_SR_res); %New RS buffer
                
                x0(Mito_free) = CaMito_basal; %initial free Ca in mitochindrium microM
                
                x0(Mito_bound1) = kmb1_on*CaMito_basal*MB1_tot/(kmb1_off+kmb1_on*CaMito_basal); %Mito buffer 1
                x0(Mito_bound2) = kmb2_on*CaMito_basal^n_mb2*MB2_tot/(kmb2_off+kmb2_on*CaMito_basal^n_mb2); %Mito buffer 2 cooperative
                
                x0(Ttubule_free) = CaTtub_basal; %µM
                x0(Ttubule_bound) = CaBTtub_basal;
                
            else
                x0(:)=x(end,:);
                
                n_p=2;
                for j=1:length(x(:,1))/n_p-1
                    x_s(j+avg,:)=mean(x((j-1)*n_p+1:(j)*n_p,:),1);
                    t_s(j+avg)=mean(t((j-1)*n_p+1:(j)*n_p));
                end
                avg=length(t_s);
                
            end
            
            if t_ini>0
                
                
                
                plot(t_s,x_s(:,Mito_free)*10,'b.',...
                    t_s,mean(x_s(:,1:m*(n-1)),2),'r')
                
                t_e=length(t_s);
                t_i=find(t_s>=t_ini-t_step+t_step/10,1);
                pt=t_s(t_i):0.0001:t_s(t_e);
                p=interp1(t_s(t_i:t_e),mean(x_s(t_i:t_e,1:m*(n-1)),2),pt);
                mid=mean(p)
                
                
                
                if rem(t_ini,3)==0
                    save last
                    pause(1)
                end
                pause(0.01)
            end
            
            t_ini
            toc
            Hz=hertz(freq);
            on=held(freq);
            optToll = odeset('RelTol', 1e-9, 'AbsTol', 1e-11); %Valori di default: 1e-3 e 1e-6, ossia 0.01% e 0.00001%
            
            
            
            [t,x]=ode15s(@(t,x) diff_diffusion_soce_auto(t,x,Hz,on,pvb_p,p_soce,anl_p,ratio_2,...
                n, m, m_SR, m_Tr, Voltofunc, CSQtofunc, MitoBtofunc, Bcitotofunc, K_D_R, K_D_L, indextofunc, mcu_syl), tspan,x0,optToll);
            
            
            
            t_ini=t_ini+t_step;
            t_end=t_end+t_step;
            
        end
        
        %clc
        t(end)
        toc
        t_mi=find(t_s>=(held(freq)+0.9-0.5),1);t_ma=find(t_s>=(held(freq)+0.9),1);
        media=mean(mean(x_s(t_mi:t_ma,1:m*(n-1)),2))
        mediasr=mean(mean(x_s(t_mi:t_ma,m*(n-1)+1:m*n),2))
        %pause
        save(filename)
        pause(0.2)
        
    end
end
