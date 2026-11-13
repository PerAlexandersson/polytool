(* Rational representatives of every returned CAD disjunct, on both sides.
   Run: wolframscript -file sample-regions.wls INPUT.json OUTPUT.json
   Representatives are retained with their full cell formula for replay. *)
ClearAll["Global`*"];
input=$ScriptCommandLine[[-2]];output=$ScriptCommandLine[[-1]];
sourceHash=IntegerString[FileHash[$InputFileName,"SHA256"],16,64];
raw=Import[input,"RawJSON"];mode=raw["mode"];
Clear[c,d];
vars=If[StringStartsQ[mode,"s2"]||StringStartsQ[mode,"circle2"],{c},
  If[MemberQ[{"s3zeroplus","s3zerominus","s3zerozero"},mode],{d},{d,c}]];
domain=Which[Length[vars]==1&&First[vars]===c,c>0,
  Length[vars]==1,d!=0,
  StringStartsQ[mode,"circle"],c>=0&&d!=0,True,c>0&&d!=0];
expr=ToExpression[raw["output"],InputForm];
sigma=Which[StringEndsQ[mode,"plus"],1,StringEndsQ[mode,"minus"],-1,True,0];
rats=Union[Flatten[Table[n/m,{n,-24,24},{m,1,24}]],Cases[expr,_Rational,Infinity],{1/36,1/27,1/54,-1/12,1/4,1/20,-1/44,-7/324,0,1,2,10,-10}];
positive=Select[rats,#>0&];
curve=Flatten[Table[
  With[{dd=(sigma*v-3)/(3*v^2),cc=(sigma*v-2)^2/v^3},
    Table[{c->cc+eps,d->dd},{eps,{0,1/100000,-1/100000}}]],{v,positive}],1];
curve=Join[curve,Flatten[Table[
  With[{dd=w^2,cc=4*sigma*w^2-8*w^3},Table[{c->cc+eps,d->dd},{eps,{0,1/100000,-1/100000}}]],
  {w,Select[positive,#<=1&]}],1]];
If[sigma==0,curve=Join[{{c->1,d->1},{c->1,d->-1},{c->0,d->1},{c->0,d->-1}},curve]];
rationalQ[x_]:=IntegerQ[x]||Head[x]===Rational;
pointValid[rules_]:=And@@(rationalQ/@(vars/.rules)) &&
  If[sigma==0&&MemberQ[vars,c],rationalQ[Sqrt[c/.rules]],True];
records={};uncovered={};
run=TimeConstrained[
 Do[
  formula=domain&&If[truth,expr,Not[expr]];
  cad=CylindricalDecomposition[formula,Evaluate[vars]];
  cells=BooleanConvert[cad,"DNF"];
  cells=If[Head[cells]===Or,List@@cells,{cells}];
  Do[
    If[cell===False,Continue[]];
    found=Quiet[TimeConstrained[FindInstance[cell,Evaluate[vars],Reals],5,{}]];
    rules=If[Length[found]>0,First[found],{}];
    If[rules==={}||!pointValid[rules]||!TrueQ[cell/.rules],
      rules=If[Length[vars]==1,
        SelectFirst[({First[vars]->#}&/@rats),TrueQ[cell/.#]&&pointValid[#]& ,{}],
        SelectFirst[curve,TrueQ[cell/.#]&&pointValid[#]&,{}]]];
    If[rules==={}&&Length[vars]==2,
      Do[
        trial=Quiet[FindInstance[(cell/.d->dd),{c},Reals]];
        If[Length[trial]>0,
          candidate=Join[First[trial],{d->dd}];
          If[pointValid[candidate],rules=candidate;Break[]]],{dd,rats}]];
    item=<|"predicate"->truth,"cell"->ToString[cell,InputForm],
      "rational"->(rules=!={}),"sample_verified"->(rules=!={}&&TrueQ[cell/.rules]),"variables"->(ToString/@vars),
      "point"->If[rules==={},{},ToString[#,InputForm]&/@(vars/.rules)]|>;
    AppendTo[records,item];If[rules==={},AppendTo[uncovered,ToString[cell,InputForm]]],{cell,cells}],
 {truth,{True,False}}],250,$Aborted];
answer=<|"schema"->"polytool.bkw-region-samples.v1","mode"->mode,"version"->$Version,
  "complete"->(run=!=$Aborted&&uncovered==={}),"samples"->records,"uncovered"->uncovered,
  "script_sha256"->sourceHash|>;
Export[output,answer,"RawJSON"];
Print[mode,": ",Length[records]," cells; uncovered=",Length[uncovered],"; complete=",answer["complete"]];
