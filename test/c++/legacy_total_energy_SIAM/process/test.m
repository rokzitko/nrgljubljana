Lambda=3;
Nnrg=1; (* Parameter Nmax in "param" of NRG Ljubljana. *)
M=Nnrg+1; (* Stevilo mest v verigi, vkljucno z 0-tim mestom. *)
PREC=20;

z=1;
NN=M-2;
Alambda = 1/2 (1+1/Lambda)/(1-1/Lambda);
factor = (1-1/Lambda)/(Log[Lambda] Sqrt[Lambda]);
Print["factor=", factor //N];
power = Lambda^(-(NN-1)/2+1-z);
omegaN = factor power;
Print["omegaN=", omegaN //N];

gammaPol = Sqrt[2 0.1/Pi Alambda];
Print["gammaPol=", gammaPol];

<<"flgreens2.m";

green[];
Print[vals];

Print[omegaN vals];

neg = Select[omegaN vals, Positive];

Print["len=", Length @ neg];
Print["sum=", 2 Plus @@ neg];

(* 
r=M/2+1;
Print["energy=", vals[[r]] ];
Print["alpha1^2=", vecs[[r, 1]]^2 ];
Print["alpha1^2/power=", vecs[[r, 1]]^2/power ];
Print["alpha1^2 Lambda^((M-3)/2)=", vecs[[r,1]]^2 Lambda^((M-3)/2) ];
*)
