syms Pl Pg Ql
display("PQ bus");
syms Irl(Vrl,Vil)
syms Iil(Vrl,Vil) 
Irl(Vrl,Vil) = (Pl*Vrl + Ql*Vil)/(Vrl^2 + Vil^2);
Iil(Vrl,Vil) = (Pl*Vil - Ql*Vrl)/(Vrl^2 + Vil^2);
dIrldVrl = diff(Irl,Vrl)
dIrldVil = diff(Irl,Vil)
dIildVrl = diff(Iil,Vrl)
dIildVil = diff(Iil,Vil)

display("PV bus");
syms Irg(Vrg,Vig,Qg)
Irg(Vrg,Vig) = (Pg*Vrg - Qg*Vig)/(Vrg^2 + Vig^2);
Iig(Vrg,Vig) = (-Pg*Vig + Qg*Vrg)/(Vrg^2 + Vig^2);
dIrgdVrg = diff(Irg,Vrg)
dIrgdVig = diff(Irg,Vig)
dIigdVrg = diff(Iig,Vrg)
dIigdVig = diff(Iig,Vig)
dIrgdQg = diff(Irg,Qg)
dIigdQg = diff(Iig,Qg)