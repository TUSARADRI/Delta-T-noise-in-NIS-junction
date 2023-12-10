# Delta-T-noise-in-NIS-junction

Subscript[\[CapitalDelta], T] noise in NIS junction for setup 1

Remove["Global`*"];

$Assumptions = 
  Z > 0 && T1 > 0 && T2 > 0 && Ef > 0 && \[CapitalDelta] > 0;
$Assumptions = 
  Ed \[Element] Reals && Z \[Element] Reals && T1 \[Element] Reals && 
   T2 \[Element] Reals;
$Assumptions = 
  Tc \[Element] Reals && 
   Ef \[Element] Reals && \[CapitalDelta] \[Element] Reals;


k = 8.617*10^-5;
Tc = 9.2;
\[CapitalDelta] = 1.764*k*Tc;

(*Fermi function for Thermovoltage*)
fe = 1/(
 1 + Exp[(Ed - V)/kT]);  fh = 1/(1 + Exp[(Ed + V)/kT]);  f1e = 1/(
 1 + Exp[(Ed - V)/kT1]); f2h = 1/(1 + Exp[(Ed + V)/kT2]);

DEfe = D[fe, Ed]; DTfe = D[fe, kT]; DEfh = D[fh, Ed]; DTfh = D[fh, kT];
DETfe = D[D[fe, kT], Ed]; DETfh = D[D[fh, kT], Ed];
DT2fe = D[fe, {kT, 2}]; DT2fh = D[fh, {kT, 2}];
DET2fe = D[D[fe, {kT, 2}], Ed]; DET2fh = D[D[fh, {kT, 2}], Ed];
DEf1e = D[f1e, Ed]; DEf2h = D[f2h, Ed];


(*Fermi function for Subscript[\[CapitalDelta], Tsh] and Subscript[\
\[CapitalDelta], Tth] noise*)

(* Without expansion of (\[CapitalDelta]T/(2 Overscript[T, _])): \
Fsh1=Simplify[(f1e-f2h)^2,Element[Ed,Reals]];
Fth1=Simplify[- (k T1)DEf1e- (k T2)DEf2h,Element[Ed,Reals]]; *)
 
Fsh = Simplify[(fe - fh)^2 + 2 (k*
\!\(\*OverscriptBox[\(T\), \(_\)]\)) (fe - fh) (DTfe + 
       DTfh) (\[CapitalDelta]T/(2 
\!\(\*OverscriptBox[\(T\), \(_\)]\))) + (k*
\!\(\*OverscriptBox[\(T\), \(_\)]\))^2 ( (DTfe - 
         DTfh)^2 + (fe - fh) (DT2fe - DT2fh)) (\[CapitalDelta]T/(2 
\!\(\*OverscriptBox[\(T\), \(_\)]\)))^2, Element[Ed, Reals]];
Fth = Simplify[- (k T1) DEfe - (k T2) DEfh + (k*
\!\(\*OverscriptBox[\(T\), \(_\)]\)) (- (k T1) DETfe + (k T2) DETfh) \
(\[CapitalDelta]T/(2 
\!\(\*OverscriptBox[\(T\), \(_\)]\))) - (k*
\!\(\*OverscriptBox[\(T\), \(_\)]\))^2/
     2 ((k T1) DET2fe + (k T2) DET2fh) (\[CapitalDelta]T/(2 
\!\(\*OverscriptBox[\(T\), \(_\)]\)))^2, Element[Ed, Reals]];


(* Equations from boundary conditions to get scattering amplitudes *)
\

a = {{0, -1, u, v}, {-1, 0, v, u}, {0, (-2*Z + I*ke), 
    I*qe*u, -I*qh*v}, {(2*Z - I*kh), 0, I*qe*v, -I*qh*u}};
b = {1, 0, I*ke + 2*Z, 0};
{ra, rb, tc, td} = 
  Refine[LinearSolve[a, b], Element[Ed, Reals]] // Simplify;

(* Coherence factors and wave vectors *)
u = (1/
   2*(1 + ((Ed^2 - \[CapitalDelta]s^2)/Ed^2)^(1/2)))^(1/2);  v = (1/
   2*(1 - ((Ed^2 - \[CapitalDelta]s^2)/Ed^2)^(1/2)))^(1/2); 

qe = (1 + (Ed^2 - \[CapitalDelta]s^2)^(1/2)/
   Ef)^(1/2); qh = (1 - (Ed^2 - \[CapitalDelta]s^2)^(1/2)/Ef)^(1/2); 
ke = (1 + Ed/Ef)^(1/2); kh = (1 - Ed/Ef)^(1/2); 

(* Reflection Probabilities *)

Reh1 = ComplexExpand[Abs[ra]^2] /. {Arg[(Ed + Ef)/Ef] -> 0, 
    Arg[1 - Ed/Ef] -> 0, Arg[1 + Ed/Ef] -> 0, Arg[1 - Ed^2/Ef^2] -> 0,
     Arg[Ed^2 - \[CapitalDelta]s^2] -> 0, 
    Arg[1 - Sqrt[Ed^2 - \[CapitalDelta]s^2]/Ef] -> 0, 
    Arg[1 + Sqrt[Ed^2 - \[CapitalDelta]s^2]/Ef] -> 0};
Ree1 = ComplexExpand[Abs[rb]^2] /. {Arg[(Ed + Ef)/Ef] -> 0, 
    Arg[1 - Ed/Ef] -> 0, Arg[1 + Ed/Ef] -> 0, Arg[1 - Ed^2/Ef^2] -> 0,
     Arg[Ed^2 - \[CapitalDelta]s^2] -> 0, 
    Arg[1 - Sqrt[Ed^2 - \[CapitalDelta]s^2]/Ef] -> 0, 
    Arg[1 + Sqrt[Ed^2 - \[CapitalDelta]s^2]/Ef] -> 0};
{Reh0, Ree0} = Abs[{ra, rb}]^2;
During evaluation of In[1]:= Remove::rmnsm: There are no symbols matching "Global`*".
In[25]:= \[CapitalDelta]
\[CapitalDelta]*Tanh[1.74*(Tc/T2 - 1)^(1/2)] /. {T2 -> 4}
Out[25]= 0.00139844
Out[26]= 0.00134652
In[27]:= Ef = 1000 \[CapitalDelta];
kT = k*
\!\(\*OverscriptBox[\(T\), \(_\)]\);
kT1 = k*T1; kT2 = k*T2;
\[CapitalDelta]T = T1 - T2;

\!\(\*OverscriptBox[\(T\), \(_\)]\) = (T1 + T2)/2;
FIi = (1 + Reh1 - Ree1);
FI = (1 + Reh0 - Ree0);
\[CapitalDelta]s = \[CapitalDelta]*Tanh[1.74*(Tc/T2 - 1)^(1/2)];
Icomp1 = Chop[Together[FIi*DTfe /. {T1 -> 6, T2 -> 4}], 10^-6];
Icomp0 = Chop[Together[FI*DTfe /. {T1 -> 6, T2 -> 4}], 10^-6];
(* Icomp2=Chop[Together[(Reh1*DTfe)/.{T1\[Rule]6,T2\[Rule]4}],10^-6];
Icomp3=Chop[Together[(Ree1*DTfe)/.{T1\[Rule]6,T2\[Rule]4}],10^-6]; *)
In[37]:= Icomp = (Icomp1 + Icomp2 - Icomp3);
In[38]:= Vthsup1 = -(2 (Ef  kT^2 (Ef^8 (16384.` + (0.` + 65536.` I) Z - 
             65536.` Z^2 + (0.` + 65536.` I) Z^3 - 
             196608.` Z^4 - (0.` + 131072.` I) Z^5) + 
          Ef^7 ((0.` - 16384.` I) + 
             98304.` Z + (0.` + 221184.` I) Z^2 - 163840.` Z^3 - 
             65536.` Z^5 - (0.` + 98304.` I) Z^6) \[CapitalDelta]s + 
          Ef^6 Z ((0.` - 36864.` I) + 
             122880.` Z + (0.` + 163840.` I) Z^2 - 
             163840.` Z^3 - (0.` + 65536.` I) Z^4 + 
             65536.` Z^5) \[CapitalDelta]s^2 + 
          Ef^5 ((0.` - 3584.` I) + 10240.` Z + (0.` + 2048.` I) Z^2 + 
             40960.` Z^3 + (0.` + 26624.` I) Z^4 - 
             32768.` Z^5) \[CapitalDelta]s^3 + 
          Ef^4 (-768.` - (0.` + 11776.` I) Z + 
             12288.` Z^2 + (0.` + 12288.` I) Z^3 - 
             4096.` Z^4) \[CapitalDelta]s^4 + 
          Ef^3 ((0.` + 384.` I) - 2304.` Z - (0.` + 3200.` I) Z^2 + 
             4096.` Z^3) \[CapitalDelta]s^5 + 
          Ef^2 (-128.` - (0.` + 512.` I) Z - 
             256.` Z^2) \[CapitalDelta]s^6 + 
          Ef ((0.` + 120.` I) - 128.` Z) \[CapitalDelta]s^7 + 
          16.` \[CapitalDelta]s^8))/(2 \[CapitalDelta]s (Ef^2 (8.` + \
(0.` + 16.` I) Z - 16.` Z^2) + 
          1.` \[CapitalDelta]s^2) (Ef^6 ((0.` - 
               1024.` I) + (2048.` + 0.` I) Z + (4096.` + 
                0.` I) Z^3 + (0.` + 4096.` I) Z^4) + 
          Ef^5 ((-512.` + 
               0.` I) - (0.` + 2048.` I) Z + (1024.` + 
                0.` I) Z^2 - (0.` + 4096.` I) Z^3 + (4096.` + 
                0.` I) Z^4 + (0.` + 4096.` I) Z^5) \[CapitalDelta]s + 
          Ef^4 Z ((-1280.` + 
               0.` I) - (0.` + 1024.` I) Z + (1024.` + 
                0.` I) Z^2 - (0.` + 
                1024.` I) Z^3) \[CapitalDelta]s^2 + 
          Ef^3 ((64.` + 
               0.` I) + (0.` + 256.` I) Z - (768.` + 
                0.` I) Z^2 - (0.` + 512.` I) Z^3) \[CapitalDelta]s^3 +
           Ef^2 ((0.` + 
               48.` I) - (64.` + 0.` I) Z + (0.` + 
                128.` I) Z^2) \[CapitalDelta]s^4 + 
          Ef ((32.` + 
               0.` I) + (0.` + 16.` I) Z) \[CapitalDelta]s^5 - (0.` + 
             4.` I) \[CapitalDelta]s^6)));
In[39]:= Dsh = (Reh1 (1 - Reh1) + Ree1 (1 - Ree1) + 2*Ree1*Reh1) Fsh;
Dth = (1 + Reh1 - Ree1) Fth;

Dsh0 = (Reh0 (1 - Reh0) + Ree0 (1 - Ree0) + 2*Ree0*Reh0) Fsh;
Dth0 = (1 + Reh0 - Ree0) Fth;
listVth = 
 ParallelTable[
  Flatten[{Z, Vthsup1 /. {T1 -> 6, T2 -> 4}}], {Z, 0.0, 3.0, 0.1}]
During evaluation of In[43]:= CloudConnect::verr: This version of Mathematica is no longer supported for use with the Wolfram Cloud beginning Sun 1 Jan 2023. Please upgrade your license or contact us.
Out[43]= {{0., 0.000385579}, {0.1, 0.000393174}, {0.2, 0.000414037}, {0.3, 
  0.000442664}, {0.4, 0.000470437}, {0.5, 0.000487878}, {0.6, 
  0.000488921}, {0.7, 0.000473901}, {0.8, 0.000447963}, {0.9, 
  0.000417089}, {1., 0.000385604}, {1.1, 0.00035586}, {1.2, 
  0.000328826}, {1.3, 0.000304705}, {1.4, 0.000283348}, {1.5, 
  0.000264469}, {1.6, 0.000247756}, {1.7, 0.00023291}, {1.8, 
  0.000219669}, {1.9, 0.000207806}, {2., 0.000197129}, {2.1, 
  0.000187478}, {2.2, 0.000178717}, {2.3, 0.000170731}, {2.4, 
  0.000163424}, {2.5, 0.000156715}, {2.6, 0.000150534}, {2.7, 
  0.000144821}, {2.8, 0.000139527}, {2.9, 0.000134607}, {3., 
  0.000130022}}
In[45]:= ListLinePlot[listVth]
Out[45]= \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[
      0.019444444444444445`], AbsoluteThickness[1.6], 
      LineBox[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQM7bbYVfXSxtJ81EwR22s/9Oici7giMf9L+s/rP
qGQla3sTYxC4bM9p+k1uLZMtVP6m/YujS/jnX7Gzh5j2wD5x/7cDVj/toeof
23vcvsARx+pgn54GAs/szSbLpe3ms4fqf2kvcOh315RoW/uzZ0Dgjf3JzqUr
p4ZaQ837YH/g7dfoaFeYez7an6w8ODEiyBxq/md7AXuD7OM9pvbnwPq/2O97
vOiOxA9jqH3f7K1vBWw5MdHInhFs3g97U9cFG/lCDaHm/bTPfFgftMzGAGre
b3sFA/1mjXY9qHv+2L9tlF91+5QO1Lx/9udlCpTP2WhD3cfgwCDD4DvxmiZU
PYODuaWc1O5JGlDzGR3m8ZydsS9H3R6sPY3JwW3RpsLjCWpQ+5gd4t12eczK
VoW6j8Uhbl6CWVWXCtQ/LA6phn5PhXYrQ81jdTjMdWdl2m8lqHvYHGZJb3m2
31MJah67g+n7XuWlSxSh5nE4KNzeFHSeS9EeAIWEzBI=

       "]]}}, {}, {}, {{}, {}}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{False, False}, {False, False}},
FrameLabel->{{None, None}, {None, None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
GridLines->{None, None},
GridLinesStyle->Directive[
GrayLevel[0.5, 0.4]],
ImagePadding->All,
Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& )}},
PlotRange->{{0, 3.0000000000000004`}, {0, 0.0004889210846446912}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.02], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\)
DshV = Dsh /. {V -> Vthsup1};
DthV = Dth /. {V -> Vthsup1};
In[48]:= listDT = ParallelTable[
   (*V=V1; Extract Vth for the current Z value*)
   
   Intsh = NIntegrate[
     DshV /. {T1 -> 6, T2 -> 4}, {Ed, -\[Infinity], \[Infinity]}, 
     Method -> {"GlobalAdaptive", "SingularityHandler" -> None, 
       "MaxErrorIncreases" -> 5000, "SymbolicProcessing" -> 0},
     MaxRecursion -> 50, WorkingPrecision -> 10]; 
   Intth = NIntegrate[
     DthV /. {T1 -> 6, T2 -> 4}, {Ed, -\[Infinity], \[Infinity]}, 
     Method -> {"GlobalAdaptive", "SingularityHandler" -> None, 
       "MaxErrorIncreases" -> 5000, "SymbolicProcessing" -> 0},
     MaxRecursion -> 50, WorkingPrecision -> 10]; 
   Inttot = Intsh + Intth; 
   DTratio = Intsh/Intth; {Z, Intsh, Intth, Inttot, 
    DTratio},(*Output Z,Intsh,Intth,Inttot,
   V*)
   {Z, 0.0, 3.0, 0.1}];
listDT = SortBy[listDT, First];(*Sort based on Z values*)
In[50]:= listDT
Out[50]= {{0., 0.00002571125765, 0.0006281700360, 0.0006538812936, 
  0.04093041084}, {0.1, 0.00002756358164, 0.0006226416725, 
  0.0006502052541, 0.04426877104}, {0.2, 0.00003399253679, 
  0.0006057145405, 0.0006397070773, 0.0561197305}, {0.3, 
  0.00004746558825, 0.0005766615425, 0.0006241271307, 
  0.0823110000}, {0.4, 0.00006988568788, 0.0005353489802, 
  0.0006052346680, 0.1305423013}, {0.5, 0.00009713582992, 
  0.0004837242247, 0.0005808600547, 0.2008082807}, {0.6, 
  0.0001174672854, 0.0004259962380, 0.0005434635234, 
  0.2757472366}, {0.7, 0.0001212186737, 0.0003673481715, 
  0.0004885668452, 0.3299830600}, {0.8, 0.0001095954574, 
  0.0003127465150, 0.0004223419724, 0.3504290285}, {0.9, 
  0.00009053162066, 0.0002637570479, 0.0003542886685, 
  0.3432386789}, {1., 0.00007097370224, 0.0002226155308, 
  0.0002935892330, 0.3188173889}, {1.1, 0.00005428206822, 
  0.0001886317049, 0.0002429137731, 0.2877674686}, {1.2, 
  0.00004122540568, 0.0001609101113, 0.0002021355170, 
  0.2562014615}, {1.3, 0.00003141746502, 0.0001383784372, 
  0.0001697959023, 0.2270401780}, {1.4, 0.00002416398566, 
  0.0001200243221, 0.0001441883077, 0.2013257417}, {1.5, 
  0.00001880982966, 0.0001049885012, 0.0001237983309, 
  0.1791608551}, {1.6, 0.00001483319177, 0.00009257793739, 
  0.0001074111292, 0.1602238307}, {1.7, 0.00001185637563, 
  0.00008224788482, 0.00009410426045, 0.1441541707}, {1.8, 
  9.597551961*10^-6, 0.00007357492766, 0.00008317247962, 
  0.1304459585}, {1.9, 7.863723378*10^-6, 0.00006623095159, 
  0.00007409467497, 0.1187318495}, {2., 6.515794833*10^-6, 
  0.00005996122002, 0.00006647701485, 0.1086668155}, {2.1, 
  5.455260189*10^-6, 0.00005456702260, 0.00006002228279, 
  0.0999735725}, {2.2, 4.611126054*10^-6, 0.00004989235371, 
  0.00005450347977, 0.0924214977}, {2.3, 3.931849486*10^-6, 
  0.00004581383269, 0.00004974568218, 0.0858223217}, {2.4, 
  3.379560505*10^-6, 0.00004223312010, 0.00004561268060, 
  0.0800215683}, {2.5, 2.926317109*10^-6, 0.00003907121442, 
  0.00004199753153, 0.0748970093}, {2.6, 2.550930950*10^-6, 
  0.00003626368062, 0.00003881461157, 0.0703439614}, {2.7, 
  2.237478309*10^-6, 0.00003375933309, 0.00003599681140, 
  0.0662773255}, {2.8, 1.973722029*10^-6, 0.00003151474264, 
  0.00003348846467, 0.0626285308}, {2.9, 1.750194129*10^-6, 
  0.00002949441041, 0.00003124460454, 0.0593398581}, {3., 
  1.559498888*10^-6, 0.00002766768002, 0.00002922717891, 
  0.0563653652}}
In[51]:= 
Export["DeltaT_setup1_NIS_Temp_64.dat", listDT];
In[52]:= DTshth1 = 
 ListLinePlot[{listDT[[All, {1, 2}]], listDT[[All, {1, 3}]], 
   listDT[[All, {1, 4}]]}, PlotRange -> {-0.00, 0.001}, 
  PlotStyle -> {{Blue, Dashed, Thickness[0.006]}, {Red, Dashed, 
     Thickness[0.006]}, {Black, Dashed, Thickness[0.006]}}, 
  AxesLabel -> {"Z", 
    "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(T\)]\)(setup 1)"}, 
  AxesStyle -> Black, 
  PlotLegends -> 
   Placed[{"\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tsh\)]\)", 
     "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tth\)]\)", 
     "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(T\)]\)"}, {Right, 
     Top}], LabelStyle -> {Bold, FontSize -> 17}]
Out[52]= \!\(
TagBox[
GraphicsBox[{{{}, {{}, {}, 
{RGBColor[0, 0, 1], PointSize[0.016666666666666666`], Thickness[
        0.006], Dashing[{Small, Small}], LineBox[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQPuX79qX/j6y27WTBDYaf8p3vkF1/M/UP5Je+kz
61anXWK0NzEGgcv2oYUzDXc+4rCHyN+0b1mz/+enQCF7iGkP7Ptu3r25pkwS
qv6x/b3Lj2dInJazT08DgWf2/0z+frlxTB6q/6W93Fy++Zd3ydifPQMCb+wF
4zPSKneLQ837YK+e2PH62ywhqPqP9nIp3SvnlvJAzf9sb73HUV1vLqv9ObD+
L/bi7pcfHK9ggNr3zV7jh19xRdhPO0aweT/sBUPvLbPc+Rnqv5/2ylsbopqk
39tBzPttH7Xqwu3wuy/sIO75Y+/9JbJ+qcITO4h5/+ybuoNenq19YAcNPofn
E3dxe4TchqpncBDZsqSh6+E1qPmMDhNrFq3VCL9sB9aexuQQ6rKk71TtBah9
zA7LhFw4eyPPQN3H4vBjl+QN854TdhD/sDh0vb8sJ5p2FGoeq0OdRs7WT0cO
Qd3D5nDoIp9mZ98BqHnsDl/0Z5xmi90LNY/D4UuFruZ/zV12AJiR3kw=
"]]}, 
{RGBColor[1, 0, 0], PointSize[0.016666666666666666`], Thickness[
        0.006], Dashing[{Small, Small}], LineBox[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQN/D8c/rJjqYj9rJgjstP9bnXZeIh3GP2n/t11f
SOqms72JMQhctt8qMNc6+akTVP6mfWWWn97VLkd7iGkP7LfL79vZvdkeqv6x
/elWzpfsr63t09NA4Jn9lL+cSSLCFlD9L+1ZUg/zstSb2J89AwJv7PcZsz+V
9jSEmvfB/tkR7Zt7dHWh6j/ap8ZfYczbqQE1/7P9Df5TIfbiqvbnwPq/2P+q
7dtpqqwEte+b/ca+n2qbyuTtGcHm/bCPuN29fUWrNNS8n/ZWKztvnXGRgJr3
2/7vnAfP+vtFoe75Y89W9mxhu6cw1Lx/9gq/J4iujxGEuo/B4frk/1du5/ND
1TM4dG/eM+/ObB6o+YwOekujLAo0uOzB2tOYHIwFBfzYWTmg9jE7bJ7C8DVP
hQ3qPhaH/Yu6tvDUsED9w+IgFSIwxYqZGWoeq4PuL6NbkpsZoe5hc+A7ft7/
WisD1Dx2B+NV5sURb//ZQczjcHDPzKz7w/TXDgDdIsru
"]]}, 
{GrayLevel[0], PointSize[0.016666666666666666`], Thickness[0.006], 
        Dashing[{Small, Small}], LineBox[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQPLi/cnqee62s+aCQI77f9MW6EW6Afjn7TfX/Kq
xuGbi72JMQhctj+9fNu5vmIXqPxNe66Etl7hq872ENMe2POelqqYzeEMVf/Y
PtF2Ifvzc4726Wkg8Mw+a+kS/1QmB6j+l/bePO1fd6y1tj97BgTe2Dd5L7T7
ZW4ONe+D/d5f5epzbI2h6j/a51qY51y9pg81/7P987pZ97rqtOzPgfV/sT/g
f3hjuqMa1L5v9rHLnilqPFOyZwSb98Peb8fCCb8tFaDm/bS/t/+mqq2GDNS8
3/Yem6IeWa2WgLrnj33/louBc86KQs37Zz/BZ/eljBxhqPsYHLiKrbONcgWh
6hkceL0uvuOt4Ieaz+hw4AufX+hkHnuw9jQmh70zn8w/IMIFtY/Z4e7VKUfZ
X7FD3cfikG/HY3SBhQ3qHxYH4cl/n86OZIGax+qwdU8na+Z9Jqh72By0tiQ8
2tfHCDWP3WHmtB/5UxIZoOZxOFxz+FMwZ+k/OwDSl8/c

         "]]}}, {}, {}, {{}, {}}, {{}, {}}}, InsetBox[
TemplateBox[{
       "\"\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \
\\(Tsh\\)]\\)\"",
        "\"\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \\(Tth\\)]\\)\
\"","\"\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \\(T\\)]\\)\""},
"LineLegend",
DisplayFunction->(FormBox[
StyleBox[
StyleBox[
PaneBox[
TagBox[
GridBox[{{
TagBox[
GridBox[{{
GraphicsBox[{{
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.054]], {
LineBox[{{0, 10}, {40, 10}}]}}, {
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.054]], {}}}, AspectRatio -> Full, ImageSize -> {40, 10}, 
                    PlotRangePadding -> None, 
                    ImagePadding -> Automatic, 
                    BaselinePosition -> (
                    Scaled[-0.11699999999999999`] -> Baseline)], #}, {
                   
GraphicsBox[{{
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
RGBColor[1, 0, 0], 
Dashing[{Small, Small}], 
Thickness[0.054]], {
LineBox[{{0, 10}, {40, 10}}]}}, {
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
RGBColor[1, 0, 0], 
Dashing[{Small, Small}], 
Thickness[0.054]], {}}}, AspectRatio -> Full, ImageSize -> {40, 10}, 
                    PlotRangePadding -> None, 
                    ImagePadding -> Automatic, 
                    BaselinePosition -> (
                    Scaled[-0.11699999999999999`] -> Baseline)], #2}, {
GraphicsBox[{{
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
GrayLevel[0], 
Dashing[{Small, Small}], 
Thickness[0.054]], {
LineBox[{{0, 10}, {40, 10}}]}}, {
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.15], 
AbsoluteThickness[1.6], 
GrayLevel[0], 
Dashing[{Small, Small}], 
Thickness[0.054]], {}}}, AspectRatio -> Full, ImageSize -> {40, 10}, 
                    PlotRangePadding -> None, 
                    ImagePadding -> Automatic, 
                    BaselinePosition -> (
                    Scaled[-0.11699999999999999`] -> Baseline)], #3}},
                   GridBoxAlignment -> {
                   "Columns" -> {Center, Left}, 
                    "Rows" -> {{Baseline}}}, AutoDelete -> False, 
                  GridBoxDividers -> {
                   "Columns" -> {{False}}, "Rows" -> {{False}}}, 
                  GridBoxItemSize -> {
                   "Columns" -> {{All}}, "Rows" -> {{All}}}, 
                  GridBoxSpacings -> {
                   "Columns" -> {{0.5}}, "Rows" -> {{0.8}}}], 
                 "Grid"]}}, 
              GridBoxAlignment -> {
               "Columns" -> {{Left}}, "Rows" -> {{Top}}}, 
              AutoDelete -> False, 
              GridBoxItemSize -> {
               "Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}, 
              
              GridBoxSpacings -> {
               "Columns" -> {{1}}, "Rows" -> {{0}}}], "Grid"], 
            Alignment -> Left, AppearanceElements -> None, 
            ImageMargins -> {{5, 5}, {5, 5}}, 
            ImageSizeAction -> "ResizeToFit"], LineIndent -> 0, 
           StripOnInput -> False], {
          FontWeight -> Bold, FontSize -> 17, FontFamily -> "Arial"}, 
          Background -> Automatic, StripOnInput -> False], 
         TraditionalForm]& ),
Editable->True,
InterpretationFunction:>(RowBox[{"LineLegend", "[", 
RowBox[{
RowBox[{"{", 
RowBox[{
RowBox[{"Directive", "[", 
RowBox[{
RowBox[{"PointSize", "[", "0.016666666666666666`", "]"}], ",", 
RowBox[{"AbsoluteThickness", "[", "1.6`", "]"}], ",", 
InterpretationBox[
ButtonBox[
TooltipBox[
GraphicsBox[{{
GrayLevel[0], 
RectangleBox[{0, 0}]}, {
GrayLevel[0], 
RectangleBox[{1, -1}]}, {
RGBColor[0, 0, 1], 
RectangleBox[{0, -1}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, 
                    FrameStyle -> RGBColor[
                    0., 0., 0.6666666666666666], FrameTicks -> None, 
                    PlotRangePadding -> None, 
                    ImageSize -> Dynamic[{
                    Automatic, 
                    1.35 CurrentValue[
                    "FontCapHeight"]/AbsoluteCurrentValue[
                    Magnification]}]], "RGBColor[0, 0, 1]"], 
                    Appearance -> None, BaseStyle -> {}, 
                    BaselinePosition -> Baseline, 
                    DefaultBaseStyle -> {}, 
                    ButtonFunction :> With[{
                    Typeset`box$ = EvaluationBox[]}, 
If[
Not[
AbsoluteCurrentValue["Deployed"]], 
                    SelectionMove[
                    Typeset`box$, All, 
                    Expression]; FrontEnd`Private`$\
ColorSelectorInitialAlpha = 1; FrontEnd`Private`$\
ColorSelectorInitialColor = RGBColor[
                    0, 0, 1]; FrontEnd`Private`$\
ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[
FrontEnd`AttachCell[Typeset`box$, 
FrontEndResource["RGBColorValueSelector"], {0, {Left, Bottom}}, {
                    Left, Top}, 
                    "ClosingActions" -> {
                    "SelectionDeparture", "ParentChanged", 
                    "EvaluatorQuit"}]]]], BaseStyle -> Inherited, 
                    Evaluator -> Automatic, Method -> "Preemptive"], 
RGBColor[0, 0, 1], Editable -> False, Selectable -> False], ",", 
RowBox[{"Dashing", "[", 
RowBox[{"{", 
RowBox[{"Small", ",", "Small"}], "}"}], "]"}], ",", 
RowBox[{"Thickness", "[", "0.006`", "]"}]}], "]"}], ",", 
RowBox[{"Directive", "[", 
RowBox[{
RowBox[{"PointSize", "[", "0.016666666666666666`", "]"}], ",", 
RowBox[{"AbsoluteThickness", "[", "1.6`", "]"}], ",", 
InterpretationBox[
ButtonBox[
TooltipBox[
GraphicsBox[{{
GrayLevel[0], 
RectangleBox[{0, 0}]}, {
GrayLevel[0], 
RectangleBox[{1, -1}]}, {
RGBColor[1, 0, 0], 
RectangleBox[{0, -1}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, 
                    FrameStyle -> RGBColor[
                    0.6666666666666666, 0., 0.], FrameTicks -> None, 
                    PlotRangePadding -> None, 
                    ImageSize -> Dynamic[{
                    Automatic, 
                    1.35 CurrentValue[
                    "FontCapHeight"]/AbsoluteCurrentValue[
                    Magnification]}]], "RGBColor[1, 0, 0]"], 
                    Appearance -> None, BaseStyle -> {}, 
                    BaselinePosition -> Baseline, 
                    DefaultBaseStyle -> {}, 
                    ButtonFunction :> With[{
                    Typeset`box$ = EvaluationBox[]}, 
If[
Not[
AbsoluteCurrentValue["Deployed"]], 
                    SelectionMove[
                    Typeset`box$, All, 
                    Expression]; FrontEnd`Private`$\
ColorSelectorInitialAlpha = 1; FrontEnd`Private`$\
ColorSelectorInitialColor = RGBColor[
                    1, 0, 0]; FrontEnd`Private`$\
ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[
FrontEnd`AttachCell[Typeset`box$, 
FrontEndResource["RGBColorValueSelector"], {0, {Left, Bottom}}, {
                    Left, Top}, 
                    "ClosingActions" -> {
                    "SelectionDeparture", "ParentChanged", 
                    "EvaluatorQuit"}]]]], BaseStyle -> Inherited, 
                    Evaluator -> Automatic, Method -> "Preemptive"], 
RGBColor[1, 0, 0], Editable -> False, Selectable -> False], ",", 
RowBox[{"Dashing", "[", 
RowBox[{"{", 
RowBox[{"Small", ",", "Small"}], "}"}], "]"}], ",", 
RowBox[{"Thickness", "[", "0.006`", "]"}]}], "]"}], ",", 
RowBox[{"Directive", "[", 
RowBox[{
RowBox[{"PointSize", "[", "0.016666666666666666`", "]"}], ",", 
RowBox[{"AbsoluteThickness", "[", "1.6`", "]"}], ",", 
InterpretationBox[
ButtonBox[
TooltipBox[
GraphicsBox[{{
GrayLevel[0], 
RectangleBox[{0, 0}]}, {
GrayLevel[0], 
RectangleBox[{1, -1}]}, {
GrayLevel[0], 
RectangleBox[{0, -1}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, 
                    FrameStyle -> GrayLevel[0.], FrameTicks -> None, 
                    PlotRangePadding -> None, 
                    ImageSize -> Dynamic[{
                    Automatic, 
                    1.35 CurrentValue[
                    "FontCapHeight"]/AbsoluteCurrentValue[
                    Magnification]}]], "GrayLevel[0]"], 
                    Appearance -> None, BaseStyle -> {}, 
                    BaselinePosition -> Baseline, 
                    DefaultBaseStyle -> {}, 
                    ButtonFunction :> With[{
                    Typeset`box$ = EvaluationBox[]}, 
If[
Not[
AbsoluteCurrentValue["Deployed"]], 
                    SelectionMove[
                    Typeset`box$, All, 
                    Expression]; FrontEnd`Private`$\
ColorSelectorInitialAlpha = 1; FrontEnd`Private`$\
ColorSelectorInitialColor = GrayLevel[
                    0]; FrontEnd`Private`$\
ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[
FrontEnd`AttachCell[Typeset`box$, 
FrontEndResource["GrayLevelColorValueSelector"], {
                    0, {Left, Bottom}}, {Left, Top}, 
                    "ClosingActions" -> {
                    "SelectionDeparture", "ParentChanged", 
                    "EvaluatorQuit"}]]]], BaseStyle -> Inherited, 
                    Evaluator -> Automatic, Method -> "Preemptive"], 
GrayLevel[0], Editable -> False, Selectable -> False], ",", 
RowBox[{"Dashing", "[", 
RowBox[{"{", 
RowBox[{"Small", ",", "Small"}], "}"}], "]"}], ",", 
RowBox[{"Thickness", "[", "0.006`", "]"}]}], "]"}]}], "}"}], ",", 
RowBox[{"{", 
RowBox[{#, ",", #2, ",", #3}], "}"}], ",", 
RowBox[{"LegendMarkers", "->", 
RowBox[{"{", 
RowBox[{
RowBox[{"{", 
RowBox[{"False", ",", "Automatic"}], "}"}], ",", 
RowBox[{"{", 
RowBox[{"False", ",", "Automatic"}], "}"}], ",", 
RowBox[{"{", 
RowBox[{"False", ",", "Automatic"}], "}"}]}], "}"}]}], ",", 
RowBox[{"Joined", "->", 
RowBox[{"{", 
RowBox[{"True", ",", "True", ",", "True"}], "}"}]}], ",", 
RowBox[{"LabelStyle", "->", 
RowBox[{"{", 
RowBox[{"Bold", ",", 
RowBox[{"FontSize", "->", "17"}]}], "}"}]}], ",", 
RowBox[{"LegendLayout", "->", "\"Column\""}]}], "]"}]& )], 
      Scaled[{0.99, 0.99}], ImageScaled[{1, 1}],
BaseStyle->{FontSize -> Larger},
FormatType->StandardForm]},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{
FormBox["\"Z\"", TraditionalForm], 
FormBox["\"\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \
\\(T\\)]\\)(setup 1)\"", TraditionalForm]},
AxesOrigin->{0, 0},
AxesStyle->GrayLevel[0],
DisplayFunction->Identity,
Frame->{{False, False}, {False, False}},
FrameLabel->{{None, None}, {None, None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
GridLines->{None, None},
GridLinesStyle->Directive[
GrayLevel[0.5, 0.4]],
ImagePadding->All,
LabelStyle->{Bold, FontSize -> 17},
Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& )}},
PlotRange->{{0, 3.0000000000000004`}, {0., 0.001}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}},
Ticks->{Automatic, Automatic}],
InterpretTemplate[Legended[
Graphics[{{}, {{{}, {}, {
Hue[0.67, 0.6, 0.6], 
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Line[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQPuX79qX/j6y27WTBDYaf8p3vkF1/M/UP5Je+kz
61anXWK0NzEGgcv2oYUzDXc+4rCHyN+0b1mz/+enQCF7iGkP7Ptu3r25pkwS
qv6x/b3Lj2dInJazT08DgWf2/0z+frlxTB6q/6W93Fy++Zd3ydifPQMCb+wF
4zPSKneLQ837YK+e2PH62ywhqPqP9nIp3SvnlvJAzf9sb73HUV1vLqv9ObD+
L/bi7pcfHK9ggNr3zV7jh19xRdhPO0aweT/sBUPvLbPc+Rnqv5/2ylsbopqk
39tBzPttH7Xqwu3wuy/sIO75Y+/9JbJ+qcITO4h5/+ybuoNenq19YAcNPofn
E3dxe4TchqpncBDZsqSh6+E1qPmMDhNrFq3VCL9sB9aexuQQ6rKk71TtBah9
zA7LhFw4eyPPQN3H4vBjl+QN854TdhD/sDh0vb8sJ5p2FGoeq0OdRs7WT0cO
Qd3D5nDoIp9mZ98BqHnsDl/0Z5xmi90LNY/D4UuFruZ/zV12AJiR3kw=
"]]}, {
Hue[0.9060679774997897, 0.6, 0.6], 
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
RGBColor[1, 0, 0], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Line[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQN/D8c/rJjqYj9rJgjstP9bnXZeIh3GP2n/t11f
SOqms72JMQhctt8qMNc6+akTVP6mfWWWn97VLkd7iGkP7LfL79vZvdkeqv6x
/elWzpfsr63t09NA4Jn9lL+cSSLCFlD9L+1ZUg/zstSb2J89AwJv7PcZsz+V
9jSEmvfB/tkR7Zt7dHWh6j/ap8ZfYczbqQE1/7P9Df5TIfbiqvbnwPq/2P+q
7dtpqqwEte+b/ca+n2qbyuTtGcHm/bCPuN29fUWrNNS8n/ZWKztvnXGRgJr3
2/7vnAfP+vtFoe75Y89W9mxhu6cw1Lx/9gq/J4iujxGEuo/B4frk/1du5/ND
1TM4dG/eM+/ObB6o+YwOekujLAo0uOzB2tOYHIwFBfzYWTmg9jE7bJ7C8DVP
hQ3qPhaH/Yu6tvDUsED9w+IgFSIwxYqZGWoeq4PuL6NbkpsZoe5hc+A7ft7/
WisD1Dx2B+NV5sURb//ZQczjcHDPzKz7w/TXDgDdIsru
"]]}, {
Hue[0.1421359549995791, 0.6, 0.6], 
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
GrayLevel[0], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Line[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQPLi/cnqee62s+aCQI77f9MW6EW6Afjn7TfX/Kq
xuGbi72JMQhctj+9fNu5vmIXqPxNe66Etl7hq872ENMe2POelqqYzeEMVf/Y
PtF2Ifvzc4726Wkg8Mw+a+kS/1QmB6j+l/bePO1fd6y1tj97BgTe2Dd5L7T7
ZW4ONe+D/d5f5epzbI2h6j/a51qY51y9pg81/7P987pZ97rqtOzPgfV/sT/g
f3hjuqMa1L5v9rHLnilqPFOyZwSb98Peb8fCCb8tFaDm/bS/t/+mqq2GDNS8
3/Yem6IeWa2WgLrnj33/louBc86KQs37Zz/BZ/eljBxhqPsYHLiKrbONcgWh
6hkceL0uvuOt4Ieaz+hw4AufX+hkHnuw9jQmh70zn8w/IMIFtY/Z4e7VKUfZ
X7FD3cfikG/HY3SBhQ3qHxYH4cl/n86OZIGax+qwdU8na+Z9Jqh72By0tiQ8
2tfHCDWP3WHmtB/5UxIZoOZxOFxz+FMwZ+k/OwDSl8/c

           "]]}}}, {}, {}, {{}, {}}, {{}, {}}}, {
      DisplayFunction -> Identity, PlotRangePadding -> {{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}}, AxesOrigin -> {0, 0}, 
       PlotRange -> {{0, 3.0000000000000004`}, {0., 0.001}}, 
       PlotRangeClipping -> True, ImagePadding -> All, 
       DisplayFunction -> Identity, AspectRatio -> GoldenRatio^(-1), 
       Axes -> {True, True}, 
       AxesLabel -> {
        "Z", "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(T\)]\)(setup \
1)"}, AxesOrigin -> {0, 0}, AxesStyle -> GrayLevel[0], 
       DisplayFunction :> Identity, 
       Frame -> {{False, False}, {False, False}}, 
       FrameLabel -> {{None, None}, {None, None}}, 
       FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}},
        GridLines -> {None, None}, GridLinesStyle -> Directive[
GrayLevel[0.5, 0.4]], LabelStyle -> {Bold, FontSize -> 17}, 
       Method -> {"CoordinatesToolOptions" -> {
          "DisplayFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& )}}, 
       PlotRange -> {{0, 3.0000000000000004`}, {0., 0.001}}, 
       PlotRangeClipping -> True, PlotRangePadding -> {{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}}, Ticks -> {Automatic, Automatic}}], 
Placed[
Unevaluated[
LineLegend[{
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
RGBColor[1, 0, 0], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Directive[
PointSize[0.016666666666666666`], 
AbsoluteThickness[1.6], 
GrayLevel[0], 
Dashing[{Small, Small}], 
Thickness[0.006]]}, {
        "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tsh\)]\)", 
         "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tth\)]\)", 
         "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(T\)]\)"}, 
        LegendMarkers -> {{False, Automatic}, {False, Automatic}, {
          False, Automatic}}, Joined -> {True, True, True}, 
        LabelStyle -> {Bold, FontSize -> 17}, 
        LegendLayout -> "Column"]], {Right, Top}, Identity]]& ],
AutoDelete->True,
Editable->True,
SelectWithContents->False,
Selectable->True]\)
In[53]:= listDT[[All, {1, 5}]]
Out[53]= {{0., 0.04093041084}, {0.1, 0.04426877104}, {0.2, 0.0561197305}, {0.3,
   0.0823110000}, {0.4, 0.1305423013}, {0.5, 0.2008082807}, {0.6, 
  0.2757472366}, {0.7, 0.3299830600}, {0.8, 0.3504290285}, {0.9, 
  0.3432386789}, {1., 0.3188173889}, {1.1, 0.2877674686}, {1.2, 
  0.2562014615}, {1.3, 0.2270401780}, {1.4, 0.2013257417}, {1.5, 
  0.1791608551}, {1.6, 0.1602238307}, {1.7, 0.1441541707}, {1.8, 
  0.1304459585}, {1.9, 0.1187318495}, {2., 0.1086668155}, {2.1, 
  0.0999735725}, {2.2, 0.0924214977}, {2.3, 0.0858223217}, {2.4, 
  0.0800215683}, {2.5, 0.0748970093}, {2.6, 0.0703439614}, {2.7, 
  0.0662773255}, {2.8, 0.0626285308}, {2.9, 0.0593398581}, {3., 
  0.0563653652}}
In[54]:= DTratio1 = 
 ListLinePlot[{listDT[[All, {1, 5}]]}, PlotRange -> {-0.00, 0.45}, 
  PlotStyle -> {{Blue, Dashed, Thickness[0.006]}, {Red, Dashed, 
     Thickness[0.006]}, {Black, Dashed, Thickness[0.006]}}, 
  AxesLabel -> {"Z", 
    "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \
\(Tsh\)]\)/\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tth\)]\)"}, 
  AxesStyle -> Black, 
  PlotLegends -> Placed[{"Setup 1"}, {Right, Top}], 
  LabelStyle -> {Bold, FontSize -> 17}]
Out[54]= \!\(
TagBox[
GraphicsBox[{{{}, {{}, {}, 
{RGBColor[0, 0, 1], PointSize[0.019444444444444445`], Thickness[
        0.006], Dashing[{Small, Small}], LineBox[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQMBZ16uv/Jlif2smSCw0/7b4d3lqauWQfkn7d2D
Z63buXuNvYkxCFy2l3zinxsqtBUqf9P+2JoiwzlbD9hDTHtgX/rD9rPolpNQ
9Y/tc1Kv7Lm+9KJ9ehoIPLM3+rjOpVDuKlT/S3u5V2YPcnOv2Z89AwJv7OXq
dGPnf78KNe+DffWfdTyN6Veg6j/ad+YutzuRdQlq/mf7co1JF2alXrA/B9b/
xT5WezHvcv6zUPu+2YvwK+/iOnrSnhFs3g/7yStP6+17dwxq3k97nea2F2ZN
R6Dm/bbnDiy4sqT4ENQ9f+wdWu8wlWw6ADXvn72c7dpHpqn7oO5jcGgoObRo
2sXdUPUMDjPW2lvcm74Taj6jg58+V977FdvtwdrTmBzsP11bUPxjK9Q+Zofn
eS+VvGu2QN3H4iBwrdi3WGcz1D8sDrbL1lzlZ9oENY/VwWjZAsbe7xug7mFz
4N9dkJzDsQFqHruDyq9z1qcT10HN43A48Y7/wuvba+wBHJXf6g==

         "]]}}, {}, {}, {{}, {}}, {{}, {}}}, InsetBox[
TemplateBox[{"\"Setup 1\""},
"LineLegend",
DisplayFunction->(FormBox[
StyleBox[
StyleBox[
PaneBox[
TagBox[
GridBox[{{
TagBox[
GridBox[{{
GraphicsBox[{{
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.175], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.054]], {
LineBox[{{0, 10}, {40, 10}}]}}, {
Directive[
EdgeForm[
Directive[
Opacity[0.3], 
GrayLevel[0]]], 
PointSize[0.175], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.054]], {}}}, AspectRatio -> Full, ImageSize -> {40, 10}, 
                    PlotRangePadding -> None, 
                    ImagePadding -> Automatic, 
                    BaselinePosition -> (
                    Scaled[-0.11699999999999999`] -> Baseline)], #}}, 
                  GridBoxAlignment -> {
                   "Columns" -> {Center, Left}, 
                    "Rows" -> {{Baseline}}}, AutoDelete -> False, 
                  GridBoxDividers -> {
                   "Columns" -> {{False}}, "Rows" -> {{False}}}, 
                  GridBoxItemSize -> {
                   "Columns" -> {{All}}, "Rows" -> {{All}}}, 
                  GridBoxSpacings -> {
                   "Columns" -> {{0.5}}, "Rows" -> {{0.8}}}], 
                 "Grid"]}}, 
              GridBoxAlignment -> {
               "Columns" -> {{Left}}, "Rows" -> {{Top}}}, 
              AutoDelete -> False, 
              GridBoxItemSize -> {
               "Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}, 
              GridBoxSpacings -> {
               "Columns" -> {{1}}, "Rows" -> {{0}}}], "Grid"], 
            Alignment -> Left, AppearanceElements -> None, 
            ImageMargins -> {{5, 5}, {5, 5}}, 
            ImageSizeAction -> "ResizeToFit"], LineIndent -> 0, 
           StripOnInput -> False], {
          FontWeight -> Bold, FontSize -> 17, FontFamily -> "Arial"}, 
          Background -> Automatic, StripOnInput -> False], 
         TraditionalForm]& ),
Editable->True,
InterpretationFunction:>(RowBox[{"LineLegend", "[", 
RowBox[{
RowBox[{"{", 
RowBox[{"Directive", "[", 
RowBox[{
RowBox[{"PointSize", "[", "0.019444444444444445`", "]"}], ",", 
RowBox[{"AbsoluteThickness", "[", "1.6`", "]"}], ",", 
InterpretationBox[
ButtonBox[
TooltipBox[
GraphicsBox[{{
GrayLevel[0], 
RectangleBox[{0, 0}]}, {
GrayLevel[0], 
RectangleBox[{1, -1}]}, {
RGBColor[0, 0, 1], 
RectangleBox[{0, -1}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, 
                    
                    FrameStyle -> RGBColor[
                    0., 0., 0.6666666666666666], FrameTicks -> None, 
                    PlotRangePadding -> None, 
                    ImageSize -> Dynamic[{
                    Automatic, 
                    1.35 CurrentValue[
                    "FontCapHeight"]/AbsoluteCurrentValue[
                    Magnification]}]], "RGBColor[0, 0, 1]"], 
                    Appearance -> None, BaseStyle -> {}, 
                    BaselinePosition -> Baseline, 
                    DefaultBaseStyle -> {}, 
                    ButtonFunction :> With[{
                    Typeset`box$ = EvaluationBox[]}, 
If[
Not[
AbsoluteCurrentValue["Deployed"]], 
                    SelectionMove[
                    Typeset`box$, All, 
                    Expression]; FrontEnd`Private`$\
ColorSelectorInitialAlpha = 1; FrontEnd`Private`$\
ColorSelectorInitialColor = RGBColor[
                    0, 0, 1]; FrontEnd`Private`$\
ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[
FrontEnd`AttachCell[Typeset`box$, 
FrontEndResource["RGBColorValueSelector"], {0, {Left, Bottom}}, {
                    Left, Top}, 
                    "ClosingActions" -> {
                    "SelectionDeparture", "ParentChanged", 
                    "EvaluatorQuit"}]]]], BaseStyle -> Inherited, 
                    Evaluator -> Automatic, Method -> "Preemptive"], 
RGBColor[0, 0, 1], Editable -> False, Selectable -> False], ",", 
RowBox[{"Dashing", "[", 
RowBox[{"{", 
RowBox[{"Small", ",", "Small"}], "}"}], "]"}], ",", 
RowBox[{"Thickness", "[", "0.006`", "]"}]}], "]"}], "}"}], ",", 
RowBox[{"{", #, "}"}], ",", 
RowBox[{"LegendMarkers", "->", 
RowBox[{"{", 
RowBox[{"{", 
RowBox[{"False", ",", "Automatic"}], "}"}], "}"}]}], ",", 
RowBox[{"Joined", "->", 
RowBox[{"{", "True", "}"}]}], ",", 
RowBox[{"LabelStyle", "->", 
RowBox[{"{", 
RowBox[{"Bold", ",", 
RowBox[{"FontSize", "->", "17"}]}], "}"}]}], ",", 
RowBox[{"LegendLayout", "->", "\"Column\""}]}], "]"}]& )], 
      Scaled[{0.99, 0.99}], ImageScaled[{1, 1}],
BaseStyle->{FontSize -> Larger},
FormatType->StandardForm]},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{
FormBox["\"Z\"", TraditionalForm], 
FormBox["\"\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \
\\(Tsh\\)]\\)/\\!\\(\\*SubscriptBox[\\(\[CapitalDelta]\\), \\(Tth\\)]\
\\)\"", TraditionalForm]},
AxesOrigin->{0, 0},
AxesStyle->GrayLevel[0],
DisplayFunction->Identity,
Frame->{{False, False}, {False, False}},
FrameLabel->{{None, None}, {None, None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
GridLines->{None, None},
GridLinesStyle->Directive[
GrayLevel[0.5, 0.4]],
ImagePadding->All,
LabelStyle->{Bold, FontSize -> 17},
Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& )}},
PlotRange->{{0, 3.0000000000000004`}, {0., 0.45}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}},
Ticks->{Automatic, Automatic}],
InterpretTemplate[Legended[
Graphics[{{}, {{{}, {}, {
Hue[0.67, 0.6, 0.6], 
Directive[
PointSize[0.019444444444444445`], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.006]], 
Line[CompressedData["
1:eJxTTMoPSmViYGCQB2IQDQMBZ16uv/Jlif2smSCw0/7b4d3lqauWQfkn7d2D
Z63buXuNvYkxCFy2l3zinxsqtBUqf9P+2JoiwzlbD9hDTHtgX/rD9rPolpNQ
9Y/tc1Kv7Lm+9KJ9ehoIPLM3+rjOpVDuKlT/S3u5V2YPcnOv2Z89AwJv7OXq
dGPnf78KNe+DffWfdTyN6Veg6j/ad+YutzuRdQlq/mf7co1JF2alXrA/B9b/
xT5WezHvcv6zUPu+2YvwK+/iOnrSnhFs3g/7yStP6+17dwxq3k97nea2F2ZN
R6Dm/bbnDiy4sqT4ENQ9f+wdWu8wlWw6ADXvn72c7dpHpqn7oO5jcGgoObRo
2sXdUPUMDjPW2lvcm74Taj6jg58+V977FdvtwdrTmBzsP11bUPxjK9Q+Zofn
eS+VvGu2QN3H4iBwrdi3WGcz1D8sDrbL1lzlZ9oENY/VwWjZAsbe7xug7mFz
4N9dkJzDsQFqHruDyq9z1qcT10HN43A48Y7/wuvba+wBHJXf6g==

           "]]}}}, {}, {}, {{}, {}}, {{}, {}}}, {
      DisplayFunction -> Identity, PlotRangePadding -> {{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}}, AxesOrigin -> {0, 0}, 
       PlotRange -> {{0, 3.0000000000000004`}, {0., 0.45}}, 
       PlotRangeClipping -> True, ImagePadding -> All, 
       DisplayFunction -> Identity, AspectRatio -> GoldenRatio^(-1), 
       Axes -> {True, True}, 
       AxesLabel -> {
        "Z", "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(Tsh\)]\)/\!\(\
\*SubscriptBox[\(\[CapitalDelta]\), \(Tth\)]\)"}, 
       AxesOrigin -> {0, 0}, AxesStyle -> GrayLevel[0], 
       DisplayFunction :> Identity, 
       Frame -> {{False, False}, {False, False}}, 
       FrameLabel -> {{None, None}, {None, None}}, 
       FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}},
        GridLines -> {None, None}, GridLinesStyle -> Directive[
GrayLevel[0.5, 0.4]], LabelStyle -> {Bold, FontSize -> 17}, 
       Method -> {"CoordinatesToolOptions" -> {
          "DisplayFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Identity[#]& )[
Part[#, 1]], 
(Identity[#]& )[
Part[#, 2]]}& )}}, 
       PlotRange -> {{0, 3.0000000000000004`}, {0., 0.45}}, 
       PlotRangeClipping -> True, PlotRangePadding -> {{
Scaled[0.02], 
Scaled[0.02]}, {0, 0}}, Ticks -> {Automatic, Automatic}}], 
Placed[
Unevaluated[
LineLegend[{
Directive[
PointSize[0.019444444444444445`], 
AbsoluteThickness[1.6], 
RGBColor[0, 0, 1], 
Dashing[{Small, Small}], 
Thickness[0.006]]}, {"Setup 1"}, 
        LegendMarkers -> {{False, Automatic}}, Joined -> {True}, 
        LabelStyle -> {Bold, FontSize -> 17}, 
        LegendLayout -> "Column"]], {Right, Top}, Identity]]& ],
AutoDelete->True,
Editable->True,
SelectWithContents->False,
Selectable->True]\)![image](https://github.com/TUSARADRI/Delta-T-noise-in-NIS-junction/assets/49065406/0cada017-d88d-404c-b87d-7c161ceef38d)
