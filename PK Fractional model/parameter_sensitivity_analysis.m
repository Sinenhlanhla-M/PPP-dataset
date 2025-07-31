(* ::Package:: *)

bigN = 500;
k12=0.6;
k21=1;
ka=2;
vm = 0.1; 
km = 0.9; 
d  = 1;  
v1 = 0.01;
F= 0.7;
c10=0;
c20=0;
c30= (F d)/v1;
ft = 10;
\[Alpha]=0.9;
\[Theta][t_]:=Exp[-k12 t]
dt=ft/bigN;
w[p_,\[Alpha]_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha];
GLSolve[parms_:{0.6,1,2,0.1,0.9,1,0.01,0.7,0.9}]:=Module[
{
c10,c20,c30,
k12,k21,ka, vm,km,d,v1,F,\[Alpha],
GLC1,
GLC2,
GLC3,
\[Omega]s
},
{k12,k21,ka, vm,km,d,v1,F,\[Alpha]}=parms;
\[Omega]s=CoefficientList[Series[Evaluate[w[1,1-\[Alpha],z]],{z,0,bigN}],z];
c10=0;
c20=0;
c30= (F d)/v1;

GLC1[0]:=GLC1[0]=c10;
GLC2[0]:=GLC2[0]=c20;
GLC3[0]:=GLC3[0]=c30;



GLC1[n_]:=GLC1[n]=GLC1[n-1]+dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])+ dt ka GLC3[n-1]-dt k12 GLC1[n-1]-dt (vm (GLC1[n-1]))/(v1 (km+GLC1[n-1]));
GLC2[n_]:=GLC2[n]=GLC2[n-1]-dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC1[n-1];
GLC3[n_]:=GLC3[n]=GLC3[n-1]-dt ka GLC3[n-1];

Return[Table[{GLC1[i],GLC2[i],GLC3[i]},{i,0,bigN}]];


]

origparms={k12,k21,ka, vm,km,d,v1,F,\[Alpha]};
\[CapitalDelta]ps=0.1origparms;
parmNum=9;
sensData=Abs[1/(2Total[UnitVector[9,parmNum]*\[CapitalDelta]ps]) (GLSolve[Evaluate[origparms+UnitVector[9,parmNum]*\[CapitalDelta]ps]]-GLSolve[Evaluate[origparms-UnitVector[9,parmNum]*\[CapitalDelta]ps]])\[Transpose]];
ts=Table[i dt,{i,0,bigN}];
sensData={ts,#}\[Transpose]&/@sensData;

doAll[parmNum_]:=(
origparms={k12,k21,ka, vm,km,d,v1,F,\[Alpha]};
\[CapitalDelta]ps=0.05origparms;
(*\[CapitalDelta]ps=0.0005origparms/origparms;*)
(*parmNum=7;*)
sensData=1/(2Total[UnitVector[9,parmNum]*\[CapitalDelta]ps]) (GLSolve[Evaluate[origparms+UnitVector[9,parmNum]*\[CapitalDelta]ps]]-GLSolve[Evaluate[origparms-UnitVector[9,parmNum]*\[CapitalDelta]ps]])\[Transpose];
ts=Table[i dt,{i,0,bigN}];
sensData={ts,#}\[Transpose]&/@sensData;
Return[ListPlot[sensData,PlotLegends->{"C1","C2","C3"},(*AxesLabel\[Rule]{"Time (Hours)","Parameter sensitivity"},*)LabelStyle->Directive[Black,Bold,18],Frame->True,FrameStyle->Directive[Bold,Black,FontSize->16],FrameLabel->{"Time (Hours)","Parameter sensitivity"},LabelStyle->Directive[Black,10],ImageSize->600,PlotRange->All,Joined->True]]
)

res=doAll/@Range[9];
res

doAll[parmNum_]:=(
origparms={k12,k21,ka, vm,km,d,v1,F,\[Alpha]};
\[CapitalDelta]ps=0.05origparms;
(*\[CapitalDelta]ps=0.0005origparms/origparms;*)
(*parmNum=7;*)
sensData=1/(2Total[UnitVector[9,parmNum]*\[CapitalDelta]ps]) (GLSolve[Evaluate[origparms+UnitVector[9,parmNum]*\[CapitalDelta]ps]]-GLSolve[Evaluate[origparms-UnitVector[9,parmNum]*\[CapitalDelta]ps]])\[Transpose];
ts=Table[i dt,{i,0,bigN}];
sensData={ts,#}\[Transpose]&/@sensData;
Return[ListPlot[sensData,PlotLegends->{"C1","C2","C3"},(*AxesLabel\[Rule]{"Time (Hours)","Parameter sensitivity"},*)LabelStyle->Directive[Black,Bold,18],Frame->True,FrameStyle->Directive[Bold,Black,FontSize->16],FrameLabel->{"Time (Hours)","Parameter sensitivity"},LabelStyle->Directive[Black,10],ImageSize->600,PlotRange->All,Joined->True]]
)

