function out = T_p(tt,T_ref,T_0,t_p)
%Fuction for the temperature at the lower side of the plate where T0 is the
%initial temperature (same as the plate) and Tref is the operating temperature 
%after the start up time tp. The monotonic function f(t) describes the increase 
%in temperature during start up.

if  (tt > t_p)
    out = T_ref;
else
    f_t = (tt/t_p).^3.*(6.*(tt/t_p).^2-15.*tt/t_p+10);
    out = (T_ref-T_0)*f_t+T_0;
end
end

