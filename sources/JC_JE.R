JC_1 <- function(Gamma_star,T_k,V_m_25,E_v_m,K_c_25,E_k_c,K_o_25,E_k_o,o_x,V_m,C_i) {
  # V_m <- V_m_25*exp(E_v_m*(T_k-298)/(8.314*298.15*T_k))
  K_c <- K_c_25*exp(E_k_c*(T_k-298)/(8.314*298.15*T_k))
  K_o <- K_o_25*exp(E_k_o*(T_k-298)/(8.314*298.15*T_k))
  O_x_over_K_o <- o_x/K_o
  Above <- C_i-Gamma_star
  Below <- C_i+K_c*(1+O_x_over_K_o)
  (V_m*(Above/Below))#*(1-exp(-k_n*LAI))/k_n
# %y=V_m.*(Above./Below).*(1-exp(-k_n*LAI)).*LAI;
}

JE_1 <- function(Gamma_star,alpha_q,I,C_i,r_JmVm,V_m_25,T_k,E_v_m) {
  # J_m <- r_JmVm*V_m_25*(airT/25) 
  J_m <- r_JmVm*V_m_25*exp(E_v_m*(T_k-298.15)/(8.314*298.15*T_k)) # %*(T_k-273)/25;
  J <- alpha_q*I*J_m/sqrt(J_m^2+alpha_q^2*I^2)
  factor_2 <- (C_i-Gamma_star)/(C_i+2*Gamma_star)
  # factor_2 <- (C_i-Gamma_star)/(4*(C_i+2*Gamma_star))
  (J/4)*factor_2#*(1-exp(-k_n*LAI))/k_n
}