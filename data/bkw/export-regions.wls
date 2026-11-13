(* Convert completed QE outputs to integer-polynomial Boolean expressions.
   Run: wolframscript -file export-regions.wls INPUT.json OUTPUT.json
   Every radical comparison is replaced with exact sign cases, then verified
   equivalent by real quantifier elimination on its physical chart. *)
ClearAll["Global`*"];
input=$ScriptCommandLine[[-2]]; output=$ScriptCommandLine[[-1]];
sourceHash=IntegerString[FileHash[$InputFileName,"SHA256"],16,64];
raw=Import[input,"RawJSON"]; mode=raw["mode"];
Clear[c,d]; vars=If[StringStartsQ[mode,"s2"]||StringStartsQ[mode,"circle2"],{c},
  If[MemberQ[{"s3zeroplus","s3zerominus","s3zerozero"},mode],{d},{c,d}]];
If[!TrueQ[raw["resolved"]],Export[output,<|"resolved"->False,"mode"->mode|>,"RawJSON"];Exit[0]];
expr=ToExpression[raw["output"],InputForm];
(* Polynomial sign is preserved by multiplying a rational numerator by its
   denominator; the denominator must separately be nonzero. *)
plain[head_,value_]:=Module[{n,den,p},
  n=Numerator[Together[value]];den=Denominator[Together[value]];p=Expand[n*den];
  If[den===1,head[p,0],den!=0 && head[p,0]]];
atom[head_,lhs_,rhs_]:=Module[{e,radicals,t,rr,f,g,ans},
  e=Together[lhs-rhs];
  radicals=DeleteDuplicates[Cases[e,Power[h_,Rational[1,2]]:>h,Infinity]];
  If[radicals==={},Return[plain[head,e]]];
  If[Length[radicals]!=1,Print["unsupported radicals: ",e];Exit[2]];
  rr=First[radicals];t=Unique["radical"];
  e=Together[e/.Power[rr,Rational[1,2]]->t];
  f=e/.t->0;g=Coefficient[e,t];
  If[!TrueQ[Expand[e-f-g*t]===0],Print["nonlinear radical: ",e];Exit[2]];
  If[head===Less,Return[atom[Greater,rhs,lhs]]];
  If[head===LessEqual,Return[atom[GreaterEqual,rhs,lhs]]];
  ans=Which[
    head===Equal,plain[Equal,f^2-g^2*rr] && plain[LessEqual,f*g],
    head===Unequal,Not[plain[Equal,f^2-g^2*rr] && plain[LessEqual,f*g]],
    head===GreaterEqual,
      (plain[GreaterEqual,g] && (plain[GreaterEqual,f] || plain[GreaterEqual,g^2*rr-f^2])) ||
      (plain[Less,g] && plain[GreaterEqual,f] && plain[GreaterEqual,f^2-g^2*rr]),
    head===Greater,
      (plain[GreaterEqual,g] && (plain[Greater,f] || plain[Greater,g^2*rr-f^2])) ||
      (plain[Less,g] && plain[Greater,f] && plain[Greater,f^2-g^2*rr]),
    True,Print["unsupported comparator: ",head];Exit[2]];
  plain[GreaterEqual,rr] && ans
];
convert[True]=True;convert[False]=False;
convert[e_And]:=And@@(convert/@List@@e);
convert[e_Or]:=Or@@(convert/@List@@e);
convert[Not[e_]]:=Not[convert[e]];
convert[e_Inequality]:=And@@Table[atom[e[[i+1]],e[[i]],e[[i+2]]],{i,1,Length[e]-2,2}];
convert[e_Equal]:=And@@Table[atom[Equal,e[[j]],e[[j+1]]],{j,Length[e]-1}];
convert[e_Unequal]:=And@@Flatten[Table[atom[Unequal,e[[i]],e[[j]]],{i,Length[e]},{j,i+1,Length[e]}]];
convert[e_]:=And@@Table[atom[Head[e],e[[j]],e[[j+1]]],{j,Length[e]-1}];
polynomial[p_]:=Module[{rules,scale},
  rules=CoefficientRules[Expand[p],vars];
  If[!And@@((IntegerQ[#]||Head[#]===Rational&)/@(Last/@rules)),Print["non-rational coefficient: ",p];Exit[4]];
  scale=Apply[LCM,Denominator[Last[#]]&/@rules];
  Map[<|"powers"->First[#],"coefficient"->ToString[scale*Last[#],InputForm]|>&,rules]
];
wire[True]=<|"op"->"true"|>;wire[False]=<|"op"->"false"|>;
wire[e_And]:=<|"op"->"and","args"->(wire/@List@@e)|>;
wire[e_Or]:=<|"op"->"or","args"->(wire/@List@@e)|>;
wire[Not[e_]]:=<|"op"->"not","arg"->wire[e]|>;
wire[e_]:=<|"op"->"sign","relation"->ToString[Head[e],InputForm],"polynomial"->polynomial[e[[1]]-e[[2]]]|>;
domain=If[MemberQ[vars,c],c>=0,True];
parts=If[Head[expr]===Or,List@@expr,{expr}];
convertedParts=convert/@parts;
converted=Or@@convertedParts;
(* Boolean closure lets us check each original disjunct independently. *)
check=TimeConstrained[And@@MapThread[
  Function[{original,translated},TrueQ[Resolve[ForAll[Evaluate[vars],
    Implies[domain,Equivalent[original,translated]]],Reals]]],
  {parts,convertedParts}],250,$Aborted];
If[!TrueQ[check],Print["equivalence check failed: ",check];Exit[3]];
stream=OpenWrite[input<>".input.txt",BinaryFormat->True];
BinaryWrite[stream,ToCharacterCode[raw["input"],"UTF8"],"Byte"];Close[stream];
answer=<|"schema"->"polytool.bkw-region-polynomials.v1","mode"->mode,
  "variables"->(ToString/@vars),"version"->$Version,
  "input_sha256"->IntegerString[FileHash[input<>".input.txt","SHA256"],16,64],
  "output_sha256"->IntegerString[FileHash[input,"SHA256"],16,64],
  "generator_sha256"->raw["script_sha256"],
  "exporter_sha256"->sourceHash,
  "equivalence_verified"->True,"formula"->wire[converted],
  "polynomial_formula"->ToString[converted,InputForm]|>;
Export[output,answer,"RawJSON"];
Print["Exported ",mode,"; ",Length[Cases[converted,_Equal|_Unequal|_Greater|_Less|_GreaterEqual|_LessEqual,Infinity]]," sign atoms"];
