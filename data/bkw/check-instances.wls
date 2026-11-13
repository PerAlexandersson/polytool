(* Independent per-instance checks of the unreduced denominator-root query.
   Run: timeout 300 wolframscript -file check-instances.wls INPUT.json OUTPUT.json
   The circle check quantifies directly over radius R and u in [-1,1].
   Every Resolve has a 25-second limit; the entire batch has a 250-second limit. *)
ClearAll["Global`*"];
Clear[u,k,h,g,a,c,d,v];

(* A tied pair has product p and sum s with s^2/p in [0,4].
   For degree three, p=u+I*v, |p|=k; Vieta determines the third root.
   Its modulus is at least the tied modulus iff c*k^3<=1.
   c is the square of the normalized second tail coefficient. *)
s2query=Exists[{u,k,h},c>0 && k>0 && k^2>u^2 && 0<=h<=4 &&
  1-2*sigma*u+sigma^2*(2*u^2-k^2)-h*c*(4*u^3-3*u*k^2)==0 &&
  -2*sigma+2*sigma^2*u-h*c*(4*u^2-k^2)==0 &&
  c*k^3<=1 && -2*u+sigma*k^2-c*k^4!=0];

(* Degree four: the tied factor is z^2-2*eta*zeta*z+zeta^2,
   zeta=u+I*v, k=|zeta|^2. Put g=b2*eta, c=b2^2, eta^2=g^2/c.
   The reciprocal quotient is monic quadratic with coefficients C,D.
   Closed-disk Schur conditions include the |D|=1 boundary separately.
   Thus the selected pair really are the two smallest roots, with multiplicity. *)
r2=2*u^2-k; r3=4*u^3-3*u*k; r4=8*u^4-8*u^2*k+k^2;
normA=k^3+d^2*k^7-2*d*k^4*r2;
normB=k^4+d^2*k^8-2*d*k^4*r4;
cross=k^3*u-2*d*k^4*r3+d^2*k^7*u;
delta=d^2*k^4;
gap=c^2*normA+4*d*c*g*cross+4*d^2*g^2*normB;
cnorm=k^3*(c^2+4*d*c*g*u+4*d^2*g^2*k);
dominance=(delta<1 && gap<=c*(1-delta)^2) ||
  (delta==1 && gap==0 && cnorm<=4*c);
s3query=Exists[{u,k,g},c>0 && k>u^2 && g>=0 && g^2<=c &&
  c*(1-sigma*r2+d*r4)-2*c*g*r3-4*d*g^2*r4==0 &&
  c*(-2*sigma*u+4*d*u*r2)-2*c*g*(4*u^2-k)-16*d*g^2*u*r2==0 &&
  g*(1+d*k*r2)+c*u*k!=0 && dominance];

(* The b2=0 chart avoids division by c. Here a=(Re zeta)^2. *)
r20=2*a-k; r40=8*a^2-8*a*k+k^2;
normB0=k^4+d^2*k^8-2*d*k^4*r40;
gap0=4*d^2*h*normB0;
s3zeroquery=Exists[{a,k,h},a>=0 && k>a && 0<h<=1 &&
  1-sigma*r20+d*(1-4*h)*r40==0 &&
  (a==0 || -sigma+2*d*(1-4*h)*r20==0) &&
  1+d*k*r20!=0 &&
  ((d^2*k^4<1 && gap0<=(1-d^2*k^4)^2) ||
   (d^2*k^4==1 && gap0==0 && 4*d^2*h*k^4<=4))];

(* Circle nonnegativity is an independent exact algebraic predicate.
   Minimize its quadratic in u on [-1,1], then use v=R^2>0. *)
endpoint=sigma*v-1+3*d*v^2;
circlequery=Exists[v,v>0 &&
  (((d<=0 || c>=16*d^2*v) && endpoint>=0 && endpoint^2>=4*c*v^3) ||
   (d>0 && c<=16*d^2*v && (4*d*sigma-c)*v-4*d-4*d^2*v^2>=0))];


input=$ScriptCommandLine[[-2]];output=$ScriptCommandLine[[-1]];
requests=Import[input,"RawJSON"];
results=TimeConstrained[Table[
  Block[{sigma=request["sigma"],cc=ToExpression[request["c"]],dd=ToExpression[request["d"]],radius,x},
    query=If[StringStartsQ[request["mode"],"circle"],
      With[{tails=ToExpression /@ request["tails"]},
        Exists[radius,radius>0 && ForAll[x,-1<=x<=1,
          Sum[tails[[j]]*radius^(j+1)*ChebyshevU[j-1,x],{j,Length[tails]}]-1>=0]]],
      If[StringStartsQ[request["mode"],"s2"],s2query,
        If[cc==0,s3zeroquery,s3query]]/.{c->cc,d->dd}];
    elapsed=AbsoluteTiming[answer=TimeConstrained[Resolve[query,Reals],25,$Aborted];][[1]];
    <|"key"->request["key"],"expected"->request["expected"],"answer"->ToString[answer,InputForm],
      "agrees"->SameQ[answer,request["expected"]],"seconds"->elapsed,
      "input"->ToString[query,InputForm],"version"->$Version|>
  ],{request,requests}],250,$Aborted];
Export[output,<|"input_sha256"->IntegerString[FileHash[input,"SHA256"],16,64],
  "script_sha256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64],"results"->results|>,"RawJSON"];
Print[results];
