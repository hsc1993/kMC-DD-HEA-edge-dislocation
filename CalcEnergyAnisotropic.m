function [W]=CalcEnergyAnisotropic(rn,links,mu,nu,a,Ec)
% this function calculates the internal energy stored in the dislocation network
% the inputs are rn, links, mu, nu, a, and Ec
% the output is energy W

TOL = 1e-9;
c11 = 1.7e5; %in MPa
c12 = 1.225e5;
c44 = 7.5e4;
bperf = (1.0/2.0).*[-1 0 1];
b = (1/6).*[-1 -1 2];
t = [-1 0 1]/norm([-1 0 1]);
bs = b*t'*t
be = cross(t,cross(b,t))
bs+be
bsnorm = norm(bs);
benorm = norm(be);
T111 = (1/sqrt(6)).*[1 -2 1; sqrt(2) sqrt(2) sqrt(2); -sqrt(3) 0 sqrt(3)];
InvT111 = inv(T111)
benorm/norm(bperf)
bsnorm/norm(bperf)

C = zeros(9,9);
C(1,1) = c11; C(2,2) = c11; C(3,3) = c11;
C(1,2) = c12; C(1,3) = c12; C(2,1) = c12; C(3,1) = c12; C(2,3) = c12; C(3,2) = c12;
C(4,4) = c44; C(4,7) = c44; C(5,5) = c44; C(5,8) = c44; C(6,6) = c44; C(6,9) = c44;
C(7,4) = c44; C(7,7) = c44; C(8,5) = c44; C(8,8) = c44; C(9,6) = c44; C(9,9) = c44;

C
Q = (1/6).*[1 2 3 -sqrt(6) -sqrt(3) sqrt(2) -sqrt(6) -sqrt(3) sqrt(2); ...
    4  2  0     0       0       -2*sqrt(2) 0   0   -2*sqrt(2); ...
    1  2  3  sqrt(6) sqrt(3) sqrt(2) sqrt(6) sqrt(3) sqrt(2); ...
    -2  2  0  sqrt(6)  0  -2*sqrt(2)  0  -2*sqrt(3)  sqrt(2); ...
    1  2  -3  -sqrt(6)  sqrt(3)  sqrt(2)  sqrt(6)  -sqrt(3)  sqrt(2); ...
    -2  2  0  0  2*sqrt(3) sqrt(2)  -sqrt(6)  0  -2*sqrt(2); ...
    -2  2  0  0  -2*sqrt(3) sqrt(2)  sqrt(6)  0  -2*sqrt(2); ...
    1  2  -3  sqrt(6)  -sqrt(3)  sqrt(2)  -sqrt(6)  sqrt(3)  sqrt(2); ...
    -2  2  0  -sqrt(6)  0  -2*sqrt(2)  0  2*sqrt(3)  sqrt(2)];

Q
Cprime = Q'*C*Q;

Cprime
for i=1:9
    for j=1:9
        if(abs(Cprime(i,j)) < TOL)
            Cprime(i,j) = 0;
        end
    end
end

H=2*c44+c12-c11;
H
A=2*c44/(c11-c12);
A
c44p = Cprime(4,4);
c55p = Cprime(5,5);
c45p = -Cprime(4,5);
Ks = sqrt(c44p*c55p-c45p^2);
Ks
mur = c44 - (1/5).*H;
mur
bsnormp = norm((1/4).*[-1 0 1]);
sigmayzp = Ks*bsnormp/(2*pi)

T110 = (1/sqrt(2)).*[0 sqrt(2) 0; -1 0 1; 1 0 1];
det(T110)
Q110 = zeros(3,3,3,3);
Qe = zeros(9,9);
Idx = [1 6 8; 9 2 4; 5 7 3];
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                Q110(k,l,i,j)=T110(i,k)*T110(j,l);
                Qe(Idx(k,l),Idx(i,j)) = Q110(k,l,i,j);
            end
        end
    end
end

Qe
CEprime = Qe'*C*Qe;
Qe'*Qe
CEprime
c11barprime = sqrt(CEprime(1,1)*CEprime(2,2));
c12prime = CEprime(1,2);
c66prime = CEprime(6,6);
c22prime = CEprime(2,2);
c11prime = CEprime(1,1);
c44pp = CEprime(4,4);
c55pp = CEprime(5,5);
M = (c11barprime+c12prime)*sqrt((c11barprime-c12prime)/(c22prime*c66prime*(c11barprime+c12prime+2*c66prime)));
M
%Kex = 
benormpp = 1.0/6.0;
sigmaxypp = -M*benormpp*c66prime/(2*pi)
Kspp = sqrt(c44pp*c55pp)

Kex = (c11barprime+c12prime)*sqrt(c66prime*(c11barprime-c12prime)/(c22prime*(c11barprime+c12prime+2*c66prime)));
Key = (c11barprime+c12prime)*sqrt(c66prime*(c11barprime-c12prime)/(c11prime*(c11barprime+c12prime+2*c66prime)));
Kex
Key
Ke = (2./3.)*Kex + Key/3
lamda_rot = c12 - (1/5).*H
mu_rot = c44 - (1/5).*H
nu_rot = lamda_rot / (2*(mu_rot+lamda_rot))
nu = c12 / (2*(c44+c11))
Iso_r = 1/(1-nu_rot);
Iso = 1/(1-nu);
Aniso_r = Ks/Ke;
Iso_r
Iso
Aniso_r
