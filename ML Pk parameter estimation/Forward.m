(* ::Package:: *)

k12=0.744668;
k21=1.5234107;
ka=2.3669796;
Vm=0.8679799;
V1=0.94791704;
Km=0.8873355;
d=0.48048388;

eq= {C1'[t] ==k21 C2[t]-k12 C1[t]+ ka C3[t]- (Vm C1[t])/(V1(Km+C1[t])),
C2'[t] == - k21 C2[t]+ k12 C1[t],

C3'[t] == -ka C3[t],
C1[0]==0,C2[0]==0,C3[0]==d

};
sol=NDSolve[eq,{C1,C2,C3},{t,0,12}];

Framed[Plot[
Evaluate[{ C1[t],C2[t],C3[t]}/.sol],
{t,0,12},
ImageSize->700,PlotRange->All,AxesLabel->{"Time (Hours)","Concentration (mg/L)"}, LabelStyle->Directive[Black,16], 
PlotLegends->{"C1-Predicted","C2-Predicted","C3-Predicted",Frame->True}
]]



