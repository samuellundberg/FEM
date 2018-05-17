function out = T_p(tt,T_ref,T_0,t_p)
if  (tt > t_p)
    out = T_ref;
else
    f_t = (tt/t_p).^3.*(6.*(tt/t_p).^2-15.*tt/t_p+10);
    out = (T_ref-T_0)*f_t+T_0;
end
end