
Tb = 290;
Tt = 220;
pb = 95000;

ent_vec = [0:0.01:10].*1e-3;
PE_vec = [0:0.01:1];
[ent_mat,PE_mat] = meshgrid(ent_vec,PE_vec);

Tb_mat = Tb.*ones(size(PE_mat));
Tt_mat = Tt.*ones(size(PE_mat));
pb_mat = pb.*ones(size(PE_mat));

[CAPE_mat,RH_mat] = calculate_CAPE_theory_const_epsilon(Tb_mat,Tt_mat,pb_mat,ent_mat,PE_mat);