
startintlab();
format long intval;
affariinit();
affariinit('ApproxMinRange');

% General parameters
global h = affari(0.5);
global t = affari(0.0);
global sim_steps = 1000;
global report_step = 20;
global iter;

% Neuron parameters
global N = affari(0.0);
global V = affari(-60.0);
global GL = affari(2.0);
global GCa = affari(4.0); % Class 1
%global GCa = affari(4.4); % Class 2
global GK = affari(8.0);
global VL = affari(-60.0);
global VCa = affari(120.0);
global VK = affari(-80.0);
global V1 = affari(-1.2);
global V2 = affari(18.0);
global V3 = affari(12.0); % Class 1
%global V3 = affari(2.0); % Class 2
global V4 = affari(17.4); % Class 1
%global V4 = affari(30.0); % Class 2
global phi = affari(1.0 / 15.0); % Class 1
%global phi = affari(1.0 / 25.0); % Class 2
global C = affari(20.0);

% Synapse parameters (excitatory)
global in_size = 50;
global R = affari(zeros(in_size, 1));
global d_R = affari(zeros(in_size, 1));
global S = affari(zeros(in_size, 1));
global I = affari(zeros(in_size, 1));
global syn_a = affari(0.25);
global syn_b = affari(0.15);
global syn_k = affari(-1.0E6);
global threshold = affari(-50.0);
global VSyn = affari(0.0);
global GSyn = affari(dlmread('conductance.dat'));
global in = dlmread('input.dat');
global in_V_lo = affari(-60.0);
global in_V_hi = affari(20.0);
global VPre = affari(zeros(in_size, 1));


% ===================== end of model parameters ======================


function d_N = dNdt ()
    global N V V3 V4 phi;

    % K+ channel activation steady-state
    % N_ss = 1 / (1 + exp(-2 (V - V3) / V4))
    temp1 = V .- V3;
    N_ss = -2 .* temp1;
    N_ss = N_ss ./ V4;
    N_ss = exp(N_ss);
    N_ss = 1 .+ N_ss;
    N_ss = 1 ./ N_ss;

    % tau of K+ channel activation
    % tau = 1 / (phi ((p + q) / 2))
    % p = exp(-(V - V3) / (2 V4))
    % q = exp( (V - V3) / (2 V4))
    temp2 = 2 .* V4;
    temp2 = temp1 / temp2;
    temp1 = -temp2;
    temp1 = exp(temp1);
    temp2 = exp(temp2);
    temp1 = temp1 .+ temp2;
    temp1 = temp1 ./ 2;
    temp1 = phi * temp1;
    temp1 = 1 ./ temp1;

    % delta of K+ channel activation
    % dN/dt = (N_ss - N) / tau
    d_N = N_ss .- N;
    d_N = d_N ./ temp1;
end


function d_V = dVdt ()
    global N S V V1 V2 VSyn VL VCa VK GSyn GL GCa GK C in_size;

    % Ca++ channel activation steady-state
    % M_ss = 1 / (1 + exp(-2 (V - V1) / V2))
    M_ss = V .- V1;
    M_ss = -2 .* M_ss;
    M_ss = M_ss ./ V2;
    M_ss = exp(M_ss);
    M_ss = 1 .+ M_ss;
    M_ss = 1 ./ M_ss;

    % Synapse current
    temp1 = VSyn - V;
    d_V = affari(0);
    for i = 1:in_size
        I(i) = temp1 * GSyn(i);
        I(i) = I(i) * S(i);
        d_V = d_V + I(i);
    end

    % Leak current
    temp1 = V .- VL;
    temp1 = temp1 .* GL;
    d_V = d_V .- temp1;

    % Ca++ current
    temp1 = V .- VCa;
    temp1 = temp1 .* GCa;
    temp1 = temp1 .* M_ss;
    d_V = d_V .- temp1;

    % K+ current
    temp1 = V .- VK;
    temp1 = temp1 .* GK;
    temp1 = temp1 .* N;
    d_V = d_V .- temp1;

    % delta of membrane potential
    % dV/dt = (I + GL (VL - V) + GCa M (VCa - V) + GK N (VK - V)) / C
    d_V = d_V ./ C;
end


function d_R = dRdt ()
    global R in in_V_lo in_V_hi VPre syn_a syn_b syn_k threshold in_size iter;

    VPre(in(iter, :) == 1) = in_V_hi;
    VPre(in(iter, :) == 0) = in_V_lo;

    % Sigmoid of threshold difference
    temp1 = VPre .- threshold;
    temp1 = temp1 .* syn_k;
    temp1 = exp(temp1);
    temp1 = temp1 .+ 1;
    temp1 = 1 ./ temp1;

    % Presynaptic transmitter release rise
    temp1 = syn_a .* temp1;

    % Presynaptic transmitter release decay
    temp2 = syn_b .* R;

    % delta of presynaptic transmitter release
    % dR/dt = a Q - b R
    % Q = 1 / (1 + e^(k(V - threshold)))
    d_R = temp1 .- temp2;
end


function d_S = dSdt ()
    global R S syn_a syn_b;

    % Postsynaptic transmitter binding rise
    temp1 = syn_a .* R;

    % Postsynaptic transmitter binding decay
    temp2 = syn_b .* S;

    % delta of postsynaptic transmitter binding
    % dS/dt = a R - b S
    d_S = temp1 .- temp2;
end


% Initialise report files
%FILE **f_time = malloc(sizeof(FILE *));;
%file_init("time", 1, f_time);

%FILE **f_nrn1_N = malloc(p_nrn1_size * sizeof(FILE *));
%file_init("nrn1_N", p_nrn1_size, f_nrn1_N);
%FILE **f_nrn1_V = malloc(p_nrn1_size * sizeof(FILE *));
%file_init("nrn1_V", p_nrn1_size, f_nrn1_V);

%FILE **f_syn_exc_R = malloc(in_size * sizeof(FILE *));
%file_init("syn_exc_R", in_size, f_syn_exc_R);
%FILE **f_syn_exc_S = malloc(in_size * sizeof(FILE *));
%file_init("syn_exc_S", in_size, f_syn_exc_S);


% Begin simulation loop
% =====================

tic;

for iter = 1:sim_steps
    if mod(iter, report_step) == 0
        disp (iter);
    end

    % Compute derivatives
    d_N = dNdt();
    d_V = dVdt();
    d_R = dRdt();
    d_S = dSdt();

    % Step system
    t = t .+ h;
    d_N = d_N .* h;
    N = N .+ d_N;
    d_V = d_V .* h;
    V = V .+ d_V;
    for j = 1:in_size
        d_R(j) = d_R(j) .* h;
        R(j) = R(j) .+ d_R(j);
        d_S(j) = d_S(j) .* h;
        S(j) = S(j) .+ d_S(j);
    end

    N
    V

    %file_write(t, 1, f_time);

    %file_write(nrn1_N, p_nrn1_size, f_nrn1_N);
    %file_write(nrn1_V, p_nrn1_size, f_nrn1_V);

    %file_write(syn_exc_R, in_size, f_syn_exc_R);
    %file_write(syn_exc_S, in_size, f_syn_exc_S);
end

toc;
