function [flux_real, flux_imag,flux_abs]=calculate_flux_density2(msh,basis,Areal,Aimag)

% calculates the flux density
%input ---data, mesh (msh), shape function info (basis),  magnetic vector potential 
% Copyright (c) 2018 Ravi Sundaria / Aalto University


flux_abs=zeros(size(basis.w,2),size(msh.t,2));
flux_real=zeros(size(basis.w,2)*2,size(msh.t,2));
flux_imag=zeros(size(basis.w,2)*2,size(msh.t,2));
if numel(Aimag)==0
    for i= 1:size(basis.w,2)
        dAnx=-1* reshape(basis.gmat(2*i-1,:),size(msh.t)).*Areal(msh.t); % delA/delnX
        By=sum(dAnx,1);
        dAny=reshape(basis.gmat(2*i,:),size(msh.t)).*Areal(msh.t); % delA/delny
        Bx=sum(dAny,1);
        flux_real(2*i-1:2*i,:)=[Bx;By];
        flux_abs(i,:)=(Bx.^2+By.^2).^(0.5);
    end
else
        for i= 1:size(basis.w,2)
        dAnx=-1* reshape(basis.gmat(2*i-1,:),size(msh.t)).*Areal(msh.t); % delA/delnX
        By=sum(dAnx,1);
        dAny=reshape(basis.gmat(2*i,:),size(msh.t)).*Areal(msh.t); % delA/delny
        Bx=sum(dAny,1);
        flux_real(2*i-1:2*i,:)=[Bx;By];
        dAnx=-1* reshape(basis.gmat(2*i-1,:),size(msh.t)).*Aimag(msh.t); % delA/delnX
        Byi=sum(dAnx,1);
        dAny=reshape(basis.gmat(2*i,:),size(msh.t)).*Aimag(msh.t); % delA/delny
        Bxi=sum(dAny,1);
        flux_imag(2*i-1:2*i,:)=[Bxi;Byi];
        flux_abs(i,:)=(Bx.^2+By.^2+Bxi.^2+Byi.^2).^(0.5);
        end
end

end