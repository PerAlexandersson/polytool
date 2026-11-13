(* Exact two-variable reduction of the quartic Schmidt-Spitzer definition.
   Run: wolframscript -file dominant-pair.wls MODE OUTPUT.json
   All tails are real; c=b2^2, d=b3 after normalization of b1=sigma.
   Use lambda=1/z, so the two smallest denominator roots become the two
   largest characteristic roots. Their factor is lambda^2-2*eta*zeta*lambda
   +zeta^2, with |eta|<=1, zeta=u+I*v and k=|zeta|^2. Vieta gives
   zeta^4-sigma*zeta^2-2*b2*eta*zeta+d*(1-4*eta^2)=0.
   Its imaginary part sets b2*eta=u*(4*u^2-2*k-sigma). Put a=u^2.
   Exact closed-unit-disk Schur inequalities check BOTH remaining roots. *)
ClearAll["Global`*"];
mode=$ScriptCommandLine[[-2]];output=$ScriptCommandLine[[-1]];
sigma=Which[StringEndsQ[mode,"plus"],1,StringEndsQ[mode,"minus"],-1,True,0];
Clear[a,k,c,d,h];
l=4*a-2*k-sigma; r2=2*a-k; r4=8*a^2-8*a*k+k^2;
equation=c*(k*(k+sigma-4*a)+d)-4*d*a*l^2;
hnorm=c^2*k+4*d*c*a*l+4*d^2*a*l^2;
gap=Expand[hnorm*(k^4+d^2)-2*d*(c^2*k^2*r2+4*d*c*a*l*k*(4*a-3*k)+4*d^2*a*l^2*r4)];
(* Reduce the Schur gap modulo Vieta before eliminating a,k. *)
schur=gap-c*(k^4-d^2)^2;
leading=k*(k+sigma-4*a)+d;
reduced=Together[Factor[PolynomialRemainder[schur,equation,c]]];
reducedSign=Expand[Numerator[reduced]*Denominator[reduced]];
(* Keep the vanishing-leading-coefficient branch explicitly. *)
dominance=(k^4>d^2 && ((leading!=0 && reducedSign<=0) ||
  (leading==0 && schur<=0))) ||
  (k^4==d^2 && gap==0 && hnorm<=4*c*k^4);
query=Exists[{a,k},c>0 && a>0 && k>a && a*l^2<=c && equation==0 &&
  l*(k^3+d*(4*a-k))+c*k!=0 && dominance];
(* With b2=0, retain eta^2=h and specialize the same closed-disk test. *)
gap0=4*d^2*h*(k^4+d^2-2*d*r4);
zeroquery=Exists[{a,k,h},a>=0 && k>a && 0<h<=1 &&
  r4-sigma*r2+d*(1-4*h)==0 && (a==0 || 2*r2==sigma) &&
  k^3+d*(4*a-k)!=0 &&
  ((k^4>d^2 && gap0<=(k^4-d^2)^2) ||
   (k^4==d^2 && gap0==0 && 4*d^2*h<=4*k^4))];
If[StringContainsQ[mode,"zero" ] && !StringEndsQ[mode,"zero"],query=zeroquery];
If[StringStartsQ[mode,"s3zero"],query=zeroquery];
elapsed=AbsoluteTiming[result=TimeConstrained[Resolve[query,Reals],250,$Aborted];][[1]];
answer=<|"schema"->"polytool.bkw-region-qe.v1","mode"->mode,"version"->$Version,
  "seconds"->elapsed,"input"->ToString[query,InputForm],
  "output"->ToString[result,InputForm],"resolved"->FreeQ[result,_Exists|_ForAll|Resolve|$Aborted],
  "script_sha256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64]|>;
Export[output,answer,"RawJSON"];
Print["Completed ",mode,"; resolved=",answer["resolved"],"; seconds=",elapsed];
