function [U,m_dot] = m_dot_function(eta,Q_in,T_HTF_in,T_HTF_out,Receiver_Height,A_Panel)

    Panel_Width=A_Panel/Receiver_Height;
    k=16.3;
    C=1522;
    Nuss=4.36;
    D_o=40/1000;
    D_th=1.25/1000;
    D_i=D_o-(D_th*2);
    

    h_tube=(Nuss*k)/D_i;

    N_tubes=floor(Panel_Width/D_o);

    R_cond=log(D_o/D_i)/(pi*Receiver_Height*k*N_tubes);
    R_conv=2/(h_tube*Receiver_Height*D_i*pi*N_tubes);

    UA=1/(R_cond+R_conv);
    U=UA/(Receiver_Height*Panel_Width)
    m_dot=(eta*Q_in)/(C*(T_HTF_out-T_HTF_in))

end

