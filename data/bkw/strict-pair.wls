(* Quartic Schmidt-Spitzer query on the actual degree-four chart c>0,d!=0.
   Run: wolframscript -file strict-pair.wls MODE OUTPUT.json
   The reduction is the same Vieta calculation as in dominant-pair.wls.
   Four characteristic roots with the same modulus sqrt(k) force
   b2=-(d/k)*Conjugate[t], hence real t. Thus an off-real tie has |d|<k^2.
   Under this strict inequality the leading coefficient used below cannot
   vanish on Vieta: its vanishing would force d=k^2. The d=0 chart is the
   separately computed cubic denominator, after trimming the trailing tail. *)
ClearAll["Global`*"];
mode=$ScriptCommandLine[[-2]];output=$ScriptCommandLine[[-1]];
sigma=Which[StringEndsQ[mode,"plus"],1,StringEndsQ[mode,"minus"],-1,True,0];
Clear[a,k,c,d];
l=4*a-2*k-sigma;r2=2*a-k;r4=8*a^2-8*a*k+k^2;
leading=k*(k+sigma-4*a)+d;
equation=c*leading-4*d*a*l^2;
hnorm=c^2*k+4*d*c*a*l+4*d^2*a*l^2;
gap=hnorm*(k^4+d^2)-2*d*(c^2*k^2*r2+4*d*c*a*l*k*(4*a-3*k)+4*d^2*a*l^2*r4)-c*(k^4-d^2)^2;
reduced=Factor[Together[PolynomialRemainder[Expand[gap],equation,c]]];
sgn=Factor[Numerator[reduced]*Denominator[reduced]];
query=Exists[{a,k},c>0 && d!=0 && a>0 && k>a && a*l^2<=c &&
  k^4>d^2 && leading!=0 && equation==0 &&
  l*(k^3+d*(4*a-k))+c*k!=0 && sgn<=0];
Print["Schur sign after Vieta: ",sgn];
elapsed=AbsoluteTiming[result=TimeConstrained[Resolve[query,Reals],250,$Aborted];][[1]];
answer=<|"schema"->"polytool.bkw-region-qe.v1","mode"->mode,"version"->$Version,
  "seconds"->elapsed,"input"->ToString[query,InputForm],
  "output"->ToString[result,InputForm],"resolved"->FreeQ[result,_Exists|_ForAll|Resolve|$Aborted],
  "script_sha256"->IntegerString[FileHash[$InputFileName,"SHA256"],16,64]|>;
Export[output,answer,"RawJSON"];
Print["Completed ",mode,"; resolved=",answer["resolved"],"; seconds=",elapsed];
