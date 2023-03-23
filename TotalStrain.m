function [ep_tot,norm_shear_strain,conv,ep_tot_conv,norm_shear_strain_conv,min]=TotalStrain(ep_inc,glideplane,ep_tot,stress_multi,norm_shear_strain_old,totaltime,total_shear,plim,j,p_file,ep_inc_conv,wp_inc_conv,ep_tot_conv,min,dt,eps_in,L,convstep)


if convstep>1000 && convstep<=5000
    eps=1.5*eps_in/(plim^3); %-6
    %disp(sprintf('increasing eps\n'));
elseif convstep>5000 && convstep<100
    eps=2*eps_in/(plim^3); %-6
    %disp(sprintf('increasing eps\n'));
else
    eps=eps_in/(plim^3); %-6
end
min_time_step=50;   %1e-12; %-12


epsilon=(ep_inc*glideplane')';
ep_normal=dot(epsilon,glideplane)*glideplane;
ep_shear=(1/2).*cross(cross(glideplane,epsilon),glideplane);

ep_tot = ep_tot + ep_inc;
ep_tot_conv_test = ep_tot_conv + ep_inc_conv;

tot_epsilon=(ep_tot*glideplane')';
tot_ep_norm=dot(tot_epsilon,glideplane)*glideplane;
shear_strain=cross(cross(glideplane,tot_epsilon),glideplane);
norm_shear_strain=norm(shear_strain);

%From last converged step
tot_epsilon_conv=(ep_tot_conv_test*glideplane')';
shear_strain_conv=cross(cross(glideplane,tot_epsilon_conv),glideplane);
norm_shear_strain_conv=norm(shear_strain_conv);
%tot_sigma=(appliedstress*glideplane')';
%tot_sigma_normal=dot(tot_sigma,glideplane)*glideplane;
%shear_stress=cross(cross(glideplane,tot_sigma),glideplane);
%norm_shear_stress=norm(shear_stress);
%norm_shear_stress_MPa=norm_shear_stress/(1.3074e-13);

%stress=0.4444444*stress_multi; %SCREWWWWW

%stress=(1/2)*stress_multi;  %EDGE

stress=stress_multi;


if(((abs(norm_shear_strain-norm_shear_strain_old)<eps)&&(convstep>min_time_step))||(norm_shear_strain>total_shear))
    conv=1;
    disp(sprintf('diff=%e eps=%e dt=%e tott=%e mint=%e norme=%e tote=%e min=%e\n',abs(norm_shear_strain-norm_shear_strain_old),eps,dt,totaltime,min_time_step,norm_shear_strain,total_shear,min));
    ep_tot_conv = ep_tot_conv + ep_inc_conv;
    fprintf(p_file,'%e\t%e\t%e\n',stress,norm_shear_strain,L);
    disp(sprintf('wrote file\n'));
else
    if (abs(norm_shear_strain-norm_shear_strain_old) < min)
        min=abs(norm_shear_strain-norm_shear_strain_old);
    end
    if(mod(j,10)==0)
        disp(sprintf('diff=%e eps=%e dt=%e tott=%e mint=%e norme=%e tote=%e min=%e step=%i\n',abs(norm_shear_strain-norm_shear_strain_old),eps,dt,totaltime,min_time_step,norm_shear_strain,total_shear,min,convstep));
    end
    conv=0;
end
