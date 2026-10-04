(* Mathematica benchmarks/nonplanarDoublePentagon.nb, d132e5365dd2a13db9cd9dbaf5c200b53d489cfd. *)
Internal={l[1],l[2]};
External={p[1],p[2],p[4],p[5]};
Propagators=#^2 & /@ {l[1],l[1]-p[5],l[1]-p[5]-p[1],l[2]-p[4]-p[2],l[2]-p[4],l[2],l[1]+l[2],l[1]+l[2]-p[1]-p[2]-p[4]-p[5],l[2]+p[5],l[1]+p[4],l[1]+p[4]+p[2]};
Replacements=Thread[{p[1]^2,p[2]^2,p[4]^2,p[5]^2,p[1] p[2],p[1] p[4],p[1] p[5],p[2] p[4],p[2] p[5],p[4] p[5]}->{0,0,0,0,s12/2,(s23-s45-s51)/2,s51/2,(-s23-s34+s51)/2,(-s12+s34-s51)/2,s45/2}];
