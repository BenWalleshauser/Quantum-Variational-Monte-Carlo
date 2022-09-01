%B.W.
%Using variational principle to estimate ground state energy

%% Parameter Selection
clear
clc
close all
%Number of protons
N = 2;

%Physical constants
e0 = 8.8541878128e-12;
hbar = 1.054571817e-34;
e_c = 1.60217663e-19;
m = 9.109e-31;
a0 = 4*pi*e0*hbar^2/(m*e_c^2);

%% Create Orthonormal Eigenfunction

%Effective nucleus charge;
Z = N-0.4;

r = sym('r',[1 N]);
theta = sym('theta',[1 N]);
phi = sym('phi',[1 N]);

%Using wavefunction assumed if no electron repulsion
psi = Z^3/(pi*a0^3)*exp(-Z*(sum(r))/a0);

%% Expectation Value of Hamiltonian

%Kinetic Energy Term
H0 = 0;
for i = 1:N
    %Laplacian
    H0 = diff(psi, r(i),2) + H0;
end 
%%
H0 = -hbar^2/(2*m)*H0;

%%
%Attractive Potential Energy Term
H1 = 0;
for i = 1:N
    H1 = 1/r(i) + H1;
end 
H1 = -N*e_c^2/(4*pi*e0)*H1*psi;

%Repulsive Potential Energy Term
H2 = 0;
for i = 1:N
    for j = 1:N
        if i == j
            continue
        end
        H2 = 1/sqrt((r(i)*cos(theta(i))-r(j)*cos(theta(j)))^2+(r(i)*sin(theta(i))*cos(phi(i))-r(j)*sin(theta(j))*cos(phi(j)))^2+(r(i)*sin(theta(i))*sin(phi(i))-r(j)*sin(theta(j))*sin(phi(j)))^2) + H2;
    end
end
H2 = H2/2;
H2 = e_c^2/(4*pi*e0)*H2*psi;

%Resulting integrand
H = H0 + H1 + H2;

%% Equilibrating Electrons

rmag = 500*a0;
del = 10;
rp = rmag.*randn(1,N);
thetap = pi*rand(1,N);
phip = 2*pi*rand(1,N);
rprime = rmag.*zeros(1,N);
thetaprime = pi*zeros(1,N);
phiprime = 2*pi*zeros(1,N);


for i = 1:N
    Set = 0;
    while Set == 0;
        rprime(i) = rp(i)+del.*a0.*randn(1,1);
        thetaprime(i) = thetap(i)+del*pi*rand(1,1);
        phiprime(i) = phip(i) + del.*2*pi*rand(1,1);
        for j = 1:N
            if i == j
                continue
            end
            rprime(j) = rp(j);
            thetaprime(j) =thetap(j);
            phiprime(j) = phip(j);
        end

        w = abs(subs(psi, [r phi theta],[rprime phiprime thetaprime])/subs(psi, [r phi theta],[rp phip thetap]))^2
        u = rand;
        v = min([w 1]);
        if u <= v
            rp = rprime;
            thetap = thetaprime;
            phip = phiprime;
            Set = 1;
        end
    end
end



%% Monte Carlo
%Choose Monte Carlo as it is a 3N dimensional integral - therefore even
%performing this integral numerically is immensely difficult

%Number of samples
NumPoints = 2000;  

%Finding average value of integrand for each Z
E_loc = zeros(NumPoints,1);

for s = 1:N:NumPoints
    for i = 1:N
        E_loc(s+i) = 6.242e18.*subs(H/psi, [r phi theta],[rp phip thetap]);
        Set = 0;
        
        %Metropolis Algorithm for calculating next configuration
        while Set == 0;
            rprime(i) = rp(i)+del.*a0.*randn(1,1);
            thetaprime(i) = thetap(i)+del*pi*rand(1,1);
            phiprime(i) = phip(i) + del.*2*pi*rand(1,1);
            for j = 1:N
                if i == j
                    continue
                end
                rprime(j) = rp(j);
                thetaprime(j) =thetap(j);
                phiprime(j) = phip(j);
            end

            w = abs(subs(psi, [r phi theta],[rprime phiprime thetaprime])/subs(psi, [r phi theta],[rp phip thetap]))^2;
            u = rand;
            v = min([w 1]);
            if u <= v
                rp = rprime;
                thetap = thetaprime;
                phip = phiprime;
                Set = 1;
            end
        end
    end
end

E = mean(E_loc);

%% Results

fprintf("\n The estimated ground state energy for N = %d is %5.3f eV",N,E);

%------ Actual Values -------
%E_helium = -79.02 eV
%E_lithium = -203.48 eV





