(* Mathematica examples/doublebox.nb, d132e5365dd2a13db9cd9dbaf5c200b53d489cfd. *)
Internal={k1,k2};
External={p1,p2,p3};
Propagators={-k1^2,-(k1+p1+p2)^2,-k2^2,-(k2+p1+p2)^2,-(k1+p1)^2,-(k1-k2)^2,-(k2-p3)^2,-(k2+p1)^2,-(k1-p3)^2};
Replacements={p1^2->0,p2^2->0,p3^2->0,p1*p2->s/2,p1*p3->t/2,p2*p3->(-s-t)/2};
