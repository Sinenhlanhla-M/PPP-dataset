(* ::Package:: *)

bigN = 500;
k12 = 0.6;
k21 = 1;
ka = 2;
vm = 0.1; 
km = 0.9; 
d  = 1;  
v1 = 0.01;
F = 0.7;
ft = 15;
\[Alpha] = 0.9;
\[Alpha]1 = 0.1;
\[Alpha]3 = 0.3;
\[Alpha]5 = 0.5;
\[Alpha]7 = 0.7;
\[Theta][t_] := Exp[-k12 t]
dt = ft/bigN;

c10 = 0;
c20 = 0;
c30 = (F d)/v1;

f[t_] := 0
g[t_] := 0
h[t_] := 0

GL
\[Alpha] = 0.1

w1[p_, \[Alpha]1_, z_] := Sum[1/j*(1 - z)^j, {j, 1, p}]^\[Alpha]1
\[Omega]s1 = CoefficientList[Series[Evaluate[w1[1, 1 - \[Alpha]1, z]], {z, 0, bigN}], z];
GLC11[0] := GLC11[0] = c10
GLC21[0] := GLC21[0] = c20
GLC31[0] := GLC31[0] = c30

GLC11[n_] := GLC11[n] = GLC11[n - 1] + dt k21^\[Alpha]1 \[Theta][n dt] (1/dt^(1 - \[Alpha]1)
Sum[\[Omega]s1[[n - j]] (GLC21[j]/\[Theta][j dt] - GLC21[0]/\[Theta][0]), {j, 0, n - 1}]) + dt ka GLC31[n - 1] - dt k12 GLC11[n - 1] - dt (vm (GLC11[n - 1]))/(v1 (km + GLC11[n - 1])) + dt f[n dt]

GLC21[n_] := GLC21[n] = GLC21[n - 1] - dt k21^\[Alpha]1 \[Theta][n dt] (1/dt^(1 - \[Alpha]1)
Sum[\[Omega]s1[[n - j]] (GLC21[j]/\[Theta][j dt] - GLC21[0]/\[Theta][0]), {j, 0, n - 1}]) + dt k12 GLC11[n - 1] +  dt g[n dt]

GLC31[n_] := GLC31[n] = GLC31[n - 1] - dt ka GLC31[n - 1] +  dt h[n dt]

GLData1 = Table[{{i dt, GLC11[i]}, {i dt, GLC21[i]}, {i dt, GLC31[i]}}, {i, 0, bigN}];
{GL11, GL21, GL31} = GLData1\[Transpose];

\[Alpha]=0.3
w3[p_,\[Alpha]3_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha]3
\[Omega]s3=CoefficientList[Series[Evaluate[w3[1,1-\[Alpha]3,z]],{z,0,bigN}],z];
GLC13[0]:=GLC13[0]=c10
GLC23[0]:=GLC23[0]=c20
GLC33[0]:=GLC33[0]=c30

GLC13[n_]:=GLC13[n]=GLC13[n-1]+dt k21^\[Alpha]3 \[Theta][n dt](1/dt^(1-\[Alpha]3) Sum[\[Omega]s3[[n-j]](GLC23[j]/\[Theta][j dt]-GLC23[0]/\[Theta][0]),{j,0,n-1}])+dt ka GLC33[n-1]-dt k12 GLC13[n-1]-dt (vm (GLC13[n-1]))/(v1 (km+GLC13[n-1])) + dt f[n dt]
GLC23[n_]:=GLC23[n]=GLC23[n-1]-dt k21^\[Alpha]3 \[Theta][n dt](1/dt^(1-\[Alpha]3) Sum[\[Omega]s3[[n-j]](GLC23[j]/\[Theta][j dt]-GLC23[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC13[n-1] +  dt g[n dt]
GLC33[n_]:=GLC33[n]=GLC33[n-1]-dt ka GLC33[n-1]+  dt h[n dt]
GLData3=Table[{{i dt,GLC13[i]},{i dt,GLC23[i]},{i dt,GLC33[i]}},{i,0,bigN}];
{GL13,GL23,GL33}=GLData3\[Transpose] ;

\[Alpha]=0.5
w5[p_,\[Alpha]5_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha]5
\[Omega]s5=CoefficientList[Series[Evaluate[w5[1,1-\[Alpha]5,z]],{z,0,bigN}],z];
GLC15[0]:=GLC15[0]=c10
GLC25[0]:=GLC25[0]=c20
GLC35[0]:=GLC53[0]=c30

GLC15[n_]:=GLC15[n]=GLC15[n-1]+dt k21^\[Alpha]5 \[Theta][n dt](1/dt^(1-\[Alpha]5) Sum[\[Omega]s5[[n-j]](GLC25[j]/\[Theta][j dt]-GLC25[0]/\[Theta][0]),{j,0,n-1}])+dt ka GLC35[n-1]-dt k12 GLC15[n-1]-dt (vm (GLC15[n-1]))/(v1 (km+GLC15[n-1])) + dt f[n dt]
GLC25[n_]:=GLC25[n]=GLC25[n-1]-dt k21^\[Alpha]5 \[Theta][n dt](1/dt^(1-\[Alpha]5) Sum[\[Omega]s5[[n-j]](GLC25[j]/\[Theta][j dt]-GLC25[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC15[n-1] +  dt g[n dt]
GLC35[n_]:=GLC35[n]=GLC35[n-1]-dt ka GLC35[n-1]+  dt h[n dt]
GLData5=Table[{{i dt,GLC15[i]},{i dt,GLC25[i]},{i dt,GLC35[i]}},{i,0,bigN}];
{GL15,GL25,GL35}=GLData5\[Transpose] ;

\[Alpha]=0.7
w7[p_,\[Alpha]7_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha]7
\[Omega]s7=CoefficientList[Series[Evaluate[w7[1,1-\[Alpha]7,z]],{z,0,bigN}],z];
GLC17[0]:=GLC17[0]=c10
GLC27[0]:=GLC27[0]=c20
GLC37[0]:=GLC37[0]=c30

GLC17[n_]:=GLC17[n]=GLC17[n-1]+dt k21^\[Alpha]7 \[Theta][n dt](1/dt^(1-\[Alpha]7) Sum[\[Omega]s7[[n-j]](GLC27[j]/\[Theta][j dt]-GLC27[0]/\[Theta][0]),{j,0,n-1}])+dt ka GLC37[n-1]-dt k12 GLC17[n-1]-dt (vm (GLC17[n-1]))/(v1 (km+GLC17[n-1])) + dt f[n dt]
GLC27[n_]:=GLC27[n]=GLC27[n-1]-dt k21^\[Alpha]7 \[Theta][n dt](1/dt^(1-\[Alpha]7) Sum[\[Omega]s7[[n-j]](GLC27[j]/\[Theta][j dt]-GLC27[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC17[n-1] +  dt g[n dt]
GLC37[n_]:=GLC37[n]=GLC37[n-1]-dt ka GLC37[n-1]+  dt h[n dt]
GLData7=Table[{{i dt,GLC17[i]},{i dt,GLC27[i]},{i dt,GLC37[i]}},{i,0,bigN}];
{GL17,GL27,GL37}=GLData7\[Transpose] ;

\[Alpha]=0.9
w[p_,\[Alpha]_,z_]:=Sum[1/j*(1-z)^j,{j,1,p}]^\[Alpha]
\[Omega]s=CoefficientList[Series[Evaluate[w[1,1-\[Alpha],z]],{z,0,bigN}],z];
GLC1[0]:=GLC1[0]=c10
GLC2[0]:=GLC2[0]=c20
GLC3[0]:=GLC3[0]=c30

GLC1[n_]:=GLC1[n]=GLC1[n-1]+dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])+dt ka GLC3[n-1]-dt k12 GLC1[n-1]-dt (vm (GLC1[n-1]))/(v1 (km+GLC1[n-1])) + dt f[n dt]
GLC2[n_]:=GLC2[n]=GLC2[n-1]-dt k21^\[Alpha] \[Theta][n dt](1/dt^(1-\[Alpha]) Sum[\[Omega]s[[n-j]](GLC2[j]/\[Theta][j dt]-GLC2[0]/\[Theta][0]),{j,0,n-1}])+dt k12 GLC1[n-1] +  dt g[n dt]
GLC3[n_]:=GLC3[n]=GLC3[n-1]-dt ka GLC3[n-1]+  dt h[n dt]
GLData=Table[{{i dt,GLC1[i]},{i dt,GLC2[i]},{i dt,GLC3[i]}},{i,0,bigN}];
{GL1,GL2,GL3}=GLData\[Transpose] ;
L1
\[Alpha]=0.1
b1[j_,\[Alpha]1_]:=dt^-\[Alpha]1/Gamma[2-\[Alpha]1] ((j+1)^(1-\[Alpha]1)-j^(1-\[Alpha]1))
L1C11[0]:=L1C11[0]=c10
L1C21[0]:=L1C21[0]=c20
L1C31[0]:=L1C31[0]=c30
L1C11[n_]:=L1C11[n]=L1C11[n-1]+dt k21^\[Alpha]1 \[Theta][n dt](Sum[b1[n-(k-1+1)-1,1-\[Alpha]1]*(L1C21[k-1+1]/\[Theta][k dt]-L1C21[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt ka L1C31[n-1]-dt k12 L1C11[n-1]-dt (vm (L1C11[n-1]))/(v1 (km+L1C11[n-1])) +  dt f[n dt];
L1C21[n_]:=L1C21[n]=L1C21[n-1]-dt k21^\[Alpha]1 \[Theta][n dt](Sum[b1[n-(k-1+1)-1,1-\[Alpha]1]*(L1C21[k-1+1]/\[Theta][k dt]-L1C21[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt k12 L1C11[n-1]+  dt g[n dt];
L1C31[n_]:=L1C31[n]=L1C31[n-1]-dt ka L1C31[n-1]+  dt h[n dt]
L1Data1=Table[{{i dt,L1C11[i]},{i dt,L1C21[i]},{i dt,L1C31[i]}},{i,0,bigN}];

\[Alpha]=0.3
b3[j_,\[Alpha]3_]:=dt^-\[Alpha]3/Gamma[2-\[Alpha]3] ((j+1)^(1-\[Alpha]3)-j^(1-\[Alpha]3))
L1C13[0]:=L1C13[0]=c10
L1C23[0]:=L1C23[0]=c20
L1C33[0]:=L1C33[0]=c30
L1C13[n_]:=L1C13[n]=L1C13[n-1]+dt k21^\[Alpha]3 \[Theta][n dt](Sum[b3[n-(k-1+1)-1,1-\[Alpha]3]*(L1C23[k-1+1]/\[Theta][k dt]-L1C23[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt ka L1C33[n-1]-dt k12 L1C13[n-1]-dt (vm (L1C13[n-1]))/(v1 (km+L1C13[n-1])) +  dt f[n dt];
L1C23[n_]:=L1C23[n]=L1C23[n-1]-dt k21^\[Alpha]3 \[Theta][n dt](Sum[b3[n-(k-1+1)-1,1-\[Alpha]3]*(L1C23[k-1+1]/\[Theta][k dt]-L1C23[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt k12 L1C13[n-1]+  dt g[n dt];
L1C33[n_]:=L1C33[n]=L1C33[n-1]-dt ka L1C33[n-1]+  dt h[n dt]
L1Data3=Table[{{i dt,L1C13[i]},{i dt,L1C23[i]},{i dt,L1C33[i]}},{i,0,bigN}];
{L113,L123,L133}=L1Data3\[Transpose] ;

\[Alpha]=0.5
b5[j_,\[Alpha]5_]:=dt^-\[Alpha]5/Gamma[2-\[Alpha]5] ((j+1)^(1-\[Alpha]5)-j^(1-\[Alpha]5))
L1C15[0]:=L1C15[0]=c10
L1C25[0]:=L1C25[0]=c20
L1C35[0]:=L1C35[0]=c30
L1C15[n_]:=L1C15[n]=L1C15[n-1]+dt k21^\[Alpha]5 \[Theta][n dt](Sum[b5[n-(k-1+1)-1,1-\[Alpha]5]*(L1C25[k-1+1]/\[Theta][k dt]-L1C25[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt ka L1C35[n-1]-dt k12 L1C15[n-1]-dt (vm (L1C15[n-1]))/(v1 (km+L1C15[n-1])) +  dt f[n dt];
L1C25[n_]:=L1C25[n]=L1C25[n-1]-dt k21^\[Alpha]5 \[Theta][n dt](Sum[b5[n-(k-1+1)-1,1-\[Alpha]5]*(L1C25[k-1+1]/\[Theta][k dt]-L1C25[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt k12 L1C15[n-1]+  dt g[n dt];
L1C35[n_]:=L1C35[n]=L1C35[n-1]-dt ka L1C35[n-1]+  dt h[n dt]
L1Data5=Table[{{i dt,L1C15[i]},{i dt,L1C25[i]},{i dt,L1C35[i]}},{i,0,bigN}];
{L115,L125,L135}=L1Data5\[Transpose] ;

\[Alpha]=0.7
b7[j_,\[Alpha]7_]:=dt^-\[Alpha]7/Gamma[2-\[Alpha]7] ((j+1)^(1-\[Alpha]7)-j^(1-\[Alpha]7))
L1C17[0]:=L1C17[0]=c10
L1C27[0]:=L1C27[0]=c20
L1C37[0]:=L1C37[0]=c30
L1C17[n_]:=L1C17[n]=L1C17[n-1]+dt k21^\[Alpha]7 \[Theta][n dt](Sum[b7[n-(k-1+1)-1,1-\[Alpha]7]*(L1C27[k-1+1]/\[Theta][k dt]-L1C27[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt ka L1C37[n-1]-dt k12 L1C17[n-1]-dt (vm (L1C17[n-1]))/(v1 (km+L1C17[n-1])) +  dt f[n dt];
L1C27[n_]:=L1C27[n]=L1C27[n-1]-dt k21^\[Alpha]7 \[Theta][n dt](Sum[b7[n-(k-1+1)-1,1-\[Alpha]7]*(L1C27[k-1+1]/\[Theta][k dt]-L1C27[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt k12 L1C17[n-1]+  dt g[n dt];
L1C37[n_]:=L1C37[n]=L1C37[n-1]-dt ka L1C37[n-1]+  dt h[n dt]
L1Data7=Table[{{i dt,L1C17[i]},{i dt,L1C27[i]},{i dt,L1C37[i]}},{i,0,bigN}];
{L117,L127,L137}=L1Data7\[Transpose] ;

\[Alpha]=0.9
b[j_,\[Alpha]_]:=dt^-\[Alpha]/Gamma[2-\[Alpha]] ((j+1)^(1-\[Alpha])-j^(1-\[Alpha]))
L1C1[0]:=L1C1[0]=c10
L1C2[0]:=L1C2[0]=c20
L1C3[0]:=L1C3[0]=c30
L1C1[n_]:=L1C1[n]=L1C1[n-1]+dt k21^\[Alpha] \[Theta][n dt](Sum[b[n-(k-1+1)-1,1-\[Alpha]]*(L1C2[k-1+1]/\[Theta][k dt]-L1C2[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt ka L1C3[n-1]-dt k12 L1C1[n-1]-dt (vm (L1C1[n-1]))/(v1 (km+L1C1[n-1])) +  dt f[n dt];
L1C2[n_]:=L1C2[n]=L1C2[n-1]-dt k21^\[Alpha] \[Theta][n dt](Sum[b[n-(k-1+1)-1,1-\[Alpha]]*(L1C2[k-1+1]/\[Theta][k dt]-L1C2[k-1]/\[Theta][(k-1) dt]),{k,1,n-1}])+dt k12 L1C1[n-1]+  dt g[n dt];
L1C3[n_]:=L1C3[n]=L1C3[n-1]-dt ka L1C3[n-1]+  dt h[n dt]

L1Data=Table[{{i dt,L1C1[i]},{i dt,L1C2[i]},{i dt,L1C3[i]}},{i,0,bigN}];
{L11,L12,L13}=L1Data\[Transpose] ;

C1 results output
ListPlot[{GL11,L111,GL13,L113,GL15,L115,GL17,L117,GL1,L11},
PlotStyle->{Blend[{Yellow,Green,Blue}],{Dashing[0.02],Red },Orange,{Dashing[0.04], Black},
{Dashing[0.03], Cyan}, {Dashing[0.01], Purple},{Dashing[0.04],Magenta},{Dashing[0.01],Green},
Blend[{Yellow,Pink,Blue}],{Dashing[0.04],Blue}},
PlotLegends->{"GL,\[Alpha]=0.1","L1,\[Alpha]=0.1","GL,\[Alpha]=0.3","L1,\[Alpha]=0.3","GL,\[Alpha]=0.5",
"L1,\[Alpha]=0.5","GL,\[Alpha]=0.7","L1,\[Alpha]=0.7","GL,\[Alpha]=0.9","L1,\[Alpha]=0.9"},
LabelStyle->Directive[Black,Bold,18],Frame->True,FrameStyle->Directive[Bold,Black,FontSize->16],
FrameLabel->{"Time (Hours)","Concentration (mg/L)"},Joined->True,ImageSize->800]

C2 results output
ListPlot[{GL21,L121,GL23,L123,GL25,L125,GL27,L127,GL2,L12},PlotStyle->{Blend[{Yellow,Green,Blue}],
{Dashing[0.02],Red },Orange,{Dashing[0.04], Black},{Dashing[0.03], Cyan}, 
{Dashing[0.01], Purple},{Dashing[0.04],Magenta},{Dashing[0.01],Green},Blend[{Yellow,Pink,Blue}],
{Dashing[0.04],Blue}},PlotLegends->{"GL,\[Alpha]=0.1","L1,\[Alpha]=0.1","GL,\[Alpha]=0.3","L1,\[Alpha]=0.3","GL,\[Alpha]=0.5",
"L1,\[Alpha]=0.5","GL,\[Alpha]=0.7","L1,\[Alpha]=0.7","GL,\[Alpha]=0.9","L1,\[Alpha]=0.9"},
LabelStyle->Directive[Black,Bold,18],Frame->True,FrameStyle->Directive[Bold,Black,FontSize->16],
FrameLabel->{"Time (Hours)","Concentration (mg/L)"},Joined->True,ImageSize->800]

C3 results output
ListPlot[{GL31,L131,GL33,L133,GL35,L135,GL37,L137,GL3,L13},
PlotStyle->{Blend[{Yellow,Green,Blue}],{Dashing[0.02],Red },Orange,{Dashing[0.04], Black},
{Dashing[0.03], Cyan}, {Dashing[0.01], Purple},
{Dashing[0.04],Magenta},{Dashing[0.01],Green},Blend[{Yellow,Pink,Blue}],
{Dashing[0.04],Blue}},
PlotLegends->{"GL,\[Alpha]=0.1","L1,\[Alpha]=0.1","GL,\[Alpha]=0.3","L1,\[Alpha]=0.3","GL,\[Alpha]=0.5","L1,\[Alpha]=0.5","GL,\[Alpha]=0.7",
"L1,\[Alpha]=0.7","GL,\[Alpha]=0.9","L1,\[Alpha]=0.9"},LabelStyle->Directive[Black,Bold,18],
Frame->True,FrameStyle->Directive[Bold,Black,FontSize->16],
FrameLabel->{"Time (Hours)","Concentration (mg/L)"},Joined->True,ImageSize->800]
