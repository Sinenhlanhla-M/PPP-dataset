(* ::Package:: *)

(* Model*)
Subscript[dC, 1]/dt=Subscript[k, 21]^-\[Alpha] \[Theta](t) D^(1-\[Alpha]) (Subscript[C, 2](t)/\[Theta](t))+Subscript[k, a] Subscript[C, 3] -Subscript[k, 12] Subscript[C, 1](t)-(Subscript[V, m]Subscript[C, 1](t))/Subscript[V, 1](Subscript[K, m]+Subscript[C, 1](t))+f(t)
Subscript[dC, 2]/dt=-Subscript[k, 21]^-\[Alpha]\[Theta](t) D^(1-\[Alpha]) (Subscript[C, 2](t)/\[Theta](t))+Subscript[k, 12] Subscript[C, 1](t)+g(t)
Subscript[dC, 3]/dt=-Subscript[k, a] Subscript[C, 3] + h(t)

Clear["Global`*"]
\[Alpha] = 9/10;
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
\[Alpha] = 9/10;
CD[x_,t_,\[Alpha]_,a_:0]:=Module[{m=Ceiling[\[Alpha]]},
Return[Assuming[{t>0,t\[Element]Reals},1/Gamma[m-\[Alpha]] Integrate[(t-\[Tau])^(m-\[Alpha]-1) Evaluate[D[x[\[Tau]],{\[Tau],m}]],{\[Tau],a,t}]]]
]
\[Theta][t_]:=Exp[-k21 t]
Clear["Global`*"]
FullSimplify/@DSolve[{
c1'[t]==k21 c2[t]-k12 c1[t]-vm c1[t]+ka c3[t],
c2'[t]==-k21 c2[t]+k12 c1[t]- vm c2[t](*+ka c3[t]*),
c3'[t]==-ka c3[t],
c1[0]==c10,
c2[0]==c20,
c3[0]==c30
},{c1[t],c2[t],c3[t]},t]
Plot[{c1[t],c2[t],c3[t]}/.%,{t,0,10}];
(*Choose exact solutions*)
C1[t_]:=1/(k12+k21) E^(-t (k12+k21)) (c20 (-1+E^((k12+k21) t)) k21+c10 (k12+E^((k12+k21) t) k21)+(c30 ka (k12 (-ka)-E^(t (k12+k21-ka)) (k12+k21) (k21-ka)+E^((k12+k21) t) k21 (k12+k21-ka)))/((ka) (k12+k21-ka)))
C2[t_]:=1/(k12+k21) E^(-t (k12+k21)) (c10 (-1+E^((k12+k21) t)) k12+c20 (E^((k12+k21) t) k12+k21)+(c30 k12 ka (-E^(t (k12+k21-ka)) (k12+k21)+ka+E^((k12+k21) t) (k12+k21-ka)))/((ka) (k12+k21-ka)))
C3[t_]:=c30 E^(-ka t)
CD[C2[#]/\[Theta][#]&,t,1-\[Alpha]]
f[t_]:=Evaluate[C1'[t]-k21^-\[Alpha] \[Theta][t]CD[C2[#]/\[Theta][#]&,t,1-\[Alpha]]+k12 C1[t]+(vm C1[t])/(v1(km+C1[t]))-ka C3[t]]
g[t_]:=Evaluate[C2'[t]+k21^-\[Alpha] \[Theta][t]CD[C2[#]/\[Theta][#]&,t,1-\[Alpha]]-k12 C1[t]]
h[t_]:=Evaluate[C3'[t]+ka C3[t]]
f[t]
g[t]
h[t]
f[t_]:=-c30 E^(t (-ka+vm)) ka+1/(k12+k21) E^(-t vm) k12 (c20 (1-E^-((k12+k21) t)) k21+c10 (E^-((k12+k21) t) k12+k21)-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 (ka-2 vm)
+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))
-(1/(k12+k21)) E^(-t vm) vm (c20 (1-E^-((k12+k21) t)) k21+c10 (E^-((k12+k21) t) k12+k21)-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 (ka-2 vm)
+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))
+(1/(k12+k21)) E^(-t vm) (c10 E^-((k12+k21) t) k12 (-k12-k21)-c20 E^-((k12+k21) t) (-k12-k21) k21-(c30 E^-((k12+k21+ka) t) (-k12-k21-ka) ka (E^(ka t) k12 (ka-2 vm)
+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm))-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 ka (ka-2 vm)
+E^(t (k12+k21+2 vm)) (k12+k21) (k12+k21+2 vm) (k21-ka+2 vm)-E^((k12+k21+ka) t) k21 (k12+k21+ka) (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))
+(E^(-t vm) vm (c20 (1-E^-((k12+k21) t)) k21+c10 (E^-((k12+k21) t) k12+k21)-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 (ka-2 vm)+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)
-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm))))/((k12+k21) v1 (km+1/(k12+k21) E^(-t vm) (c20 (1-E^-((k12+k21) t)) k21+c10 (E^-((k12+k21) t) k12+k21)
-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 (ka-2 vm)+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))))
-(1/(k21^(9/10) (k12+k21) Gamma[9/10])) E^(-k21 t) ((c10 E^(t (k21-vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)+(c20 E^(t (k21-vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)
-(c10 E^(t (k21-vm)) k12 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)-(c20 E^(t (k21-vm)) k12 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)
+(c30 E^(t (k21-vm)) k12^2 k21 ka (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 k21^2 ka (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-vm)) k12 k21 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 ka^2 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))
-(2 c30 E^(t (k21-vm)) k12 ka vm^2 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 k21 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((k21-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((k21-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))+(c10 E^(-t (k12+vm)) k12^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)
-(c20 E^(-t (k12+vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)+(c10 E^(-t (k12+vm)) k12 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)
-(c20 E^(-t (k12+vm)) k21 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)+(c30 E^(-t (k12+vm)) k12^2 ka^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
-(2 c30 E^(-t (k12+vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))+(c30 E^(-t (k12+vm)) k12 ka^2 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
-(2 c30 E^(-t (k12+vm)) k12 ka vm^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))-(c30 E^(t (k21-ka+vm)) k12^2 k21 ka (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
-(c30 E^(t (k21-ka+vm)) k12 k21^2 ka (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-ka+vm)) k12^2 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-ka+vm)) k12 k21 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))-(c30 E^(t (k21-ka+vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
-(c30 E^(t (k21-ka+vm)) k12 k21 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm)))

g[t_]:=-(1/(k12+k21)) E^(-t vm) vm (c10 (1-E^-((k12+k21) t)) k12+c20 (k12+E^-((k12+k21) t) k21)+(c30 E^-((k12+k21+ka) t) k12 ka (-E^(t (k12+k21+2 vm)) (k12+k21)+E^(ka t) (ka-2 vm)+E^((k12+k21+ka) t) (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))
-(1/(k12+k21)) E^(-t vm) k12 (c20 (1-E^-((k12+k21) t)) k21+c10 (E^-((k12+k21) t) k12+k21)-(c30 E^-((k12+k21+ka) t) ka (E^(ka t) k12 (ka-2 vm)+E^(t (k12+k21+2 vm)) (k12+k21) (k21-ka+2 vm)
-E^((k12+k21+ka) t) k21 (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))+1/(k12+k21) E^(-t vm) (-c10 E^-((k12+k21) t) k12 (-k12-k21)+c20 E^-((k12+k21) t) (-k12-k21) k21
+(c30 E^-((k12+k21+ka) t) k12 (-k12-k21-ka) ka (-E^(t (k12+k21+2 vm)) (k12+k21)+E^(ka t) (ka-2 vm)+E^((k12+k21+ka) t) (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm))+
(c30 E^-((k12+k21+ka) t) k12 ka (E^(ka t) ka (ka-2 vm)-E^(t (k12+k21+2 vm)) (k12+k21) (k12+k21+2 vm)+E^((k12+k21+ka) t) (k12+k21+ka) (k12+k21-ka+2 vm)))/((ka-2 vm) (k12+k21-ka+2 vm)))
+(1/(k21^(9/10) (k12+k21) Gamma[9/10])) E^(-k21 t) ((c10 E^(t (k21-vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)+(c20 E^(t (k21-vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)
-(c10 E^(t (k21-vm)) k12 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)-(c20 E^(t (k21-vm)) k12 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/(k21-vm)^(9/10)
+(c30 E^(t (k21-vm)) k12^2 k21 ka (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 k21^2 ka (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-vm)) k12 k21 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 ka^2 vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))
-(2 c30 E^(t (k21-vm)) k12 ka vm^2 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((ka-2 vm) (k21-vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-vm)) k12 k21 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((k21-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-vm)]))/((k21-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))+(c10 E^(-t (k12+vm)) k12^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)
-(c20 E^(-t (k12+vm)) k12 k21 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)+(c10 E^(-t (k12+vm)) k12 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)
-(c20 E^(-t (k12+vm)) k21 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/(-k12-vm)^(9/10)+(c30 E^(-t (k12+vm)) k12^2 ka^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
-(2 c30 E^(-t (k12+vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))+(c30 E^(-t (k12+vm)) k12 ka^2 vm (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))
-(2 c30 E^(-t (k12+vm)) k12 ka vm^2 (Gamma[9/10]-Gamma[9/10,-t (k12+vm)]))/((-k12-vm)^(9/10) (-ka+2 vm) (k12+k21-ka+2 vm))-(c30 E^(t (k21-ka+vm)) k12^2 k21 ka (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
-(c30 E^(t (k21-ka+vm)) k12 k21^2 ka (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))+(c30 E^(t (k21-ka+vm)) k12^2 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
+(c30 E^(t (k21-ka+vm)) k12 k21 ka^2 (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))-(c30 E^(t (k21-ka+vm)) k12^2 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm))
-(c30 E^(t (k21-ka+vm)) k12 k21 ka vm (Gamma[9/10]-Gamma[9/10,t (k21-ka+vm)]))/((ka-2 vm) (k21-ka+vm)^(9/10) (k12+k21-ka+2 vm)))

h[t_]:=c30 E^(t (-ka+vm)) ka+c30 E^(t (-ka+vm)) (-ka+vm)
FullSimplify[f[t]]
FullSimplify[g[t]]
FullSimplify[h[t]]

