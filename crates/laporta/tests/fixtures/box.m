(* Literal FIRE family definition; squared propagators may have either sign. *)
Internal = {k};
External = {p1,p2,p3};
Propagators = {-k^2, -(k+p1)^2, -(k+p1+p2)^2, -(k+p1+p2+p3)^2};
Replacements = {p1^2->0,p2^2->0,p3^2->0,p1*p2->s/2,p2*p3->t/2,p1*p3->(-s-t)/2};
