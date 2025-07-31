(* ::Package:: *)

(*Fractional model*)
ndsol=NDSolve[{
C1'[t]==k21 CaputoD[C2[t],{t,1-\[Alpha]}]+ka C3[t]-k12 C1[t]-(vm C1[t])/(v1(km+C1[t])),
C2'[t]==-k21 CaputoD[C2[t],{t,1-\[Alpha]}]+k12 C1[t],
C3'[t]==-ka C3[t],
C1[0]==c10,
C2[0]==c20,
C3[0]==c30
},{C1[t],C2[t],C3[t]},{t,0,ft}]

ndsol=NDSolve[{C1'[t]==k21 C2[t]+ka C3[t]-k12 C1[t]-(vm C1[t])/(v1(km+C1[t])),
C2'[t]==-k21 C2[t]+k12 C1[t],
C3'[t]==-ka C3[t],
C1[0]==c10,
C2[0]==c20,
C3[0]==c30
},{C1[t],C2[t],C3[t]},{t,0,ft}]

bigN = 500;
getData[datasheet_]:=(Module[{dataset,res},
labels=datasheet[[1]];
dataset=Dataset[Association/@(Thread[datasheet[[1,All]]->#]&/@datasheet[[2;;-1]])];
tids=Normal[DeleteDuplicates[dataset[All,"tid"]]];
tids=DeleteCases[tids,_?(#==""&)];
times=Normal[DeleteDuplicates[dataset[All,"p_time"]]];
times=DeleteCases[times,_?(#==""&)];
(*res=Table[Normal[dataset[Select[#tid==tid&],labels[[-1]]]],{tid,tids}];*)
res=Table[Normal[dataset[Select[(#day==day&&#tid==tid&)],labels[[-1]]]],{tid,tids}];
Return[res]
]
)
(*Data Prelim*)
data=Import["C:\\P1\\p1\\Simulation\\Oral dose\\Data fittting\\CuratedData.xlsx"];
labels=data[[2,1]];
dataset=Dataset[Association/@(Thread[data[[2]][[1,All]]->#]&/@data[[2]][[2;;-1]])];
tids=Normal[DeleteDuplicates[dataset[All,"tid"]]];
tids=DeleteCases[tids,_?(#==""&)];
res=Table[Normal[dataset[Select[#tid==tid&&#day==3&],labels[[-1]]]],{tid,tids}];
times=Normal[DeleteDuplicates[dataset[All,"p_time"]]];
times=DeleteCases[times,_?(#==""&)];

(*Data Fitting*)
getData[datasheet_]:=(Module[{dataset,res},
labels=datasheet[[1]];
dataset=Dataset[Association/@(Thread[datasheet[[1,All]]->#]&/@datasheet[[2;;-1]])];
tids=Normal[DeleteDuplicates[dataset[All,"tid"]]];
tids=DeleteCases[tids,_?(#==""&)];
times=Normal[DeleteDuplicates[dataset[All,"p_time"]]];
times=DeleteCases[times,_?(#==""&)];
(*res=Table[Normal[dataset[Select[#tid==tid&],labels[[-1]]]],{tid,tids}];*)
res=Table[Normal[dataset[Select[(#day==day&&#tid==tid&)],labels[[-1]]]],{tid,tids}];
Return[res]
]
)


getSol[parms_]:=Module[(*{k12, k21,ka,vm,km,F,v1,\[Alpha],c10,c20,c30}*){},
{k12, k21,ka,vm,km,F,v1,\[Alpha]}=parms;
c10=Mean[res][[1]];(*predose value*)
c20=0;
c30=F*d/v1;

dt=ft/bigN;
pos=Flatten[Position[N@Range[0,ft,dt],_?(MemberQ[times,#]&)]];

w[p_,\[Alpha]_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha];
\[Omega]s=CoefficientList[Series[Evaluate[w[1,1-\[Alpha],z]],{z,0,bigN}],z];
GLC1[0]:=GLC1[0]=c10;
GLC2[0]:=GLC2[0]=c20;
GLC3[0]:=GLC3[0]=c30;

GLC1[n_]:=GLC1[n]=GLC1[n-1]+dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])-dt k12 GLC1[n-1]+dt ka c30 Exp[-ka n dt]-dt (vm (GLC1[n-1]))/(v1 (km+GLC1[n-1])) ;
GLC2[n_]:=GLC2[n]=GLC2[n-1]-dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC1[n-1];
GLC3[n_]:=GLC3[n]=GLC3[n-1]-dt ka c30 Exp[-ka n dt];

GLData=Table[{{i dt,GLC1[i]},{i dt,GLC2[i]},{i dt,GLC3[i]}},{i,0,bigN}];
{GL1,GL2,GL3}=GLData\[Transpose] ;

Return[GL1[[pos]]];
];

getObj[parms_]:=(
ClearAll[GLC1,GLC2,GLC3];
Return[Norm[getSol[parms]-{times,Mean[res]}\[Transpose],1]]
)

data=Import["C:\\P1\\p1\\Simulation\\Oral dose\\Data fittting\\CuratedData.xlsx"];
drugData="LPV";
drugData="RTV";
day=1;
If[drugData=="LPV",
res=getData[data[[2]]]; (*LPV data*)
res=DeleteCases[res,_?(#=={}&)];
If[day==1,d=0.4];
If[day==2,d=0.4]; 
If[day==3,d=0.6];
If[day==4,d=0.8];
,
res=getData[data[[4]]]; (*RTV data*)
res=DeleteCases[res,_?(#=={}&)];
If[day==1,d=0.1];
If[day==2,d=0.1]; 
If[day==3,d=0.15];
If[day==4,d=0.2];
]

Print["Using data from: "<>drugData<>"\nDay: "<>ToString[day]<>"\nDose is: "<>ToString[d]]
ft=12;
bigN=120;
\[Theta][t_]:=Exp[-k12 t]
f[t_]:=0
g[t_]:=0

{k12, k21,ka,vm,km,F,v1,\[Alpha]}={0.6,1,2,0.1,0.9,0.7,0.01,0.5};
parms={k12, k21,ka,vm,km,F,v1,\[Alpha]};
best=getObj[parms]
p=0.1;
count=0;

(*Parameter algorithm*)
Monitor[
While[best>=0.5&&p>0.04,
oldbest=best;
For[j=1,j<=Length[parms],j++,
candp=getObj[parms*(ConstantArray[1,Length[parms]]+p*UnitVector[Length[parms],j])];
candn=getObj[parms*(ConstantArray[1,Length[parms]]-p*UnitVector[Length[parms],j])];
If[candp<best||candn<best,
If[candp<candn,
parms=parms*(ConstantArray[1,Length[parms]]+p*UnitVector[Length[parms],j]);
best=candp;
,
parms=parms*(ConstantArray[1,Length[parms]]-p*UnitVector[Length[parms],j]);
best=candn;
];
];
];
If[oldbest==best,count=count+1,count=0];
If[count>1,p=0.9*p;count=0;];
If[RandomReal[]<0.1,p=p*1.2;];

];
,{ListLinePlot[{getSol[parms],{times,Mean[res]}\[Transpose]}],best,p,parms}]
parms
best
ListLinePlot[{getSol[parms],{times,Mean[res]}\[Transpose]},AxesLabel->{"Time (Hours)","Concentration (mg/L)"},PlotLegends->{"Mathematical Model","Clinical Data"},
PlotMarkers->{Automatic,14},
ImageSize->800,LabelStyle->Directive[Black,10]]


