function Efv_work = Efv_work(MU,NU,rn,b,a,appliedstress,omega_tensor)

mu=MU;
nu=NU;

x1= [0 0 0];
x2= [0 0 100];

x1= rn(1:end-1,1:3);
x2= rn(2:end,1:3);
x = 0.5*(x1+x2);

% Stress(i,:)=[s_11 s_22 s_33 s_12 s_23 s_13]
[fieldPointStress] = FieldPointStress(x,x1,x2,b,a,mu,nu); % same as [mu] -> N/A^2
numTensor = size(fieldPointStress,1);
stressTensor = zeros(3,3,numTensor);
work = zeros(numTensor,1);
work_ev = zeros(numTensor,1);

for ii = 1:numTensor
    tensorTemp = [fieldPointStress(ii,1) fieldPointStress(ii,4) fieldPointStress(ii,6);
                fieldPointStress(ii,4) fieldPointStress(ii,2) fieldPointStress(ii,5);
                fieldPointStress(ii,6) fieldPointStress(ii,5) fieldPointStress(ii,3)];
    stressTensor(:,:,ii) = tensorTemp+appliedstress

    for kk = 1:3
        for jj = 1:3
            work(ii) = work(ii)-stressTensor(kk,jj,ii)*omega_tensor(kk,jj); % N*A   
            % stresstensor>0 -> region tensile -> work>0 needed for vacancy
            % formation -> omega_tensor_v<0 and minus sign is added.
        end
    end
    work_ev(ii) = work(ii)*1e-10*6.242e+18 %ev
end
Efv_work = max(work_ev);

end