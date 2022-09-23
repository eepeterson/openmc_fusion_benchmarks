iter_port_plug_fwcadis
c
c Attila calculation "ITER_Port_Plug_FWCADIS"
c no motivating tally so we get even errors everywhere
c -------------------------- Input Information ------------------------------ 80
c Attila GUI created MCNP6 Input
c Attila Version 10.2.1.5984
c Input File Creation Date: Fri Sep 23 15:35:03 2022
c RxMesher version          : 1.0.0
c Simmetrix MeshSim version : 14.0-200321
c Global mesh size          : 0.5 m
c Generated                 : 2022-09-23T15:12:24-04:00
c Solid Geometry :
c   Filename     : ITER_port_plug_benchmark.x_t
c   Last changed : 2022-09-23T15:00:28-04:00
c   MD5 checksum : cf915e6a86b6e05eeb36343e1c7e07b2
c Note: RTT Mesh has added cell flags for MCNP Abaqus part and pseudo-cell.
c Associated Abaqus Unstructured Mesh :
c   Filename     : iter_port_plug_benchmark.abaq
c   Generated    : 2022-09-23T15:17:48
c   MD5 checksum : cd0652d512738cf8117f51bf21068162
c n_points = 15706
c n_sides  = 14833
c n_cells  = 88547
c Mesh Bounding Box (cm):
c  x: -150.000 -  150.000
c  y: -150.000 -  150.000
c  z:    0.000 -  700.000
c
c Number of Attila Regions                 : 15
c Number of Abaqus Parts/MCNP Pseudo-Cells : 15
c Number of Materials                      : 2
c
c  Mesh Region/Pseudo-Cell Information
c   Attila Region #    : 1
c   Attila Region Name : "Source"
c   Abaqus Part #      : 1
c   Abaqus Part Name   : "Source"
c   MCNP Pseudo-cell # : 1
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 312955 cm**3
c     # Cells          : 269
c     % of Total Cells : 0.30% 
c
c   Attila Region #    : 2
c   Attila Region Name : "Frame"
c   Abaqus Part #      : 2
c   Abaqus Part Name   : "Frame"
c   MCNP Pseudo-cell # : 2
c   Material           : "Steel-el_1"
c    MCNP Material     : m1
c    Density           : 8 g/cc
c   Mesh Data
c     Meshed Volume    : 1.29183e+07 cm**3
c     # Cells          : 23642
c     % of Total Cells : 26.70% 
c
c   Attila Region #    : 3
c   Attila Region Name : "Plug"
c   Abaqus Part #      : 3
c   Abaqus Part Name   : "Plug"
c   MCNP Pseudo-cell # : 3
c   Material           : "Steel-el_1"
c    MCNP Material     : m1
c    Density           : 8 g/cc
c   Mesh Data
c     Meshed Volume    : 1.4779e+06 cm**3
c     # Cells          : 18500
c     % of Total Cells : 20.89% 
c
c   Attila Region #    : 4
c   Attila Region Name : "BackPlate"
c   Abaqus Part #      : 4
c   Abaqus Part Name   : "BackPlate"
c   MCNP Pseudo-cell # : 4
c   Material           : "Steel-el_1"
c    MCNP Material     : m1
c    Density           : 8 g/cc
c   Mesh Data
c     Meshed Volume    : 108207 cm**3
c     # Cells          : 248
c     % of Total Cells : 0.28% 
c
c   Attila Region #    : 5
c   Attila Region Name : "Void0"
c   Abaqus Part #      : 5
c   Abaqus Part Name   : "Void0"
c   MCNP Pseudo-cell # : 5
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 3.1315e+06 cm**3
c     # Cells          : 3308
c     % of Total Cells : 3.74% 
c
c   Attila Region #    : 6
c   Attila Region Name : "Void1"
c   Abaqus Part #      : 6
c   Abaqus Part Name   : "Void1"
c   MCNP Pseudo-cell # : 6
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 36995.4 cm**3
c     # Cells          : 9372
c     % of Total Cells : 10.58% 
c
c   Attila Region #    : 7
c   Attila Region Name : "Void2"
c   Abaqus Part #      : 7
c   Abaqus Part Name   : "Void2"
c   MCNP Pseudo-cell # : 7
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 337790 cm**3
c     # Cells          : 8548
c     % of Total Cells : 9.65% 
c
c   Attila Region #    : 8
c   Attila Region Name : "Void3"
c   Abaqus Part #      : 8
c   Abaqus Part Name   : "Void3"
c   MCNP Pseudo-cell # : 8
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 2.34458e+06 cm**3
c     # Cells          : 7828
c     % of Total Cells : 8.84% 
c
c   Attila Region #    : 9
c   Attila Region Name : "Void4"
c   Abaqus Part #      : 9
c   Abaqus Part Name   : "Void4"
c   MCNP Pseudo-cell # : 9
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 940278 cm**3
c     # Cells          : 2723
c     % of Total Cells : 3.08% 
c
c   Attila Region #    : 10
c   Attila Region Name : "Tally0"
c   Abaqus Part #      : 10
c   Abaqus Part Name   : "Tally0"
c   MCNP Pseudo-cell # : 10
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 7044.87 cm**3
c     # Cells          : 308
c     % of Total Cells : 0.35% 
c
c   Attila Region #    : 11
c   Attila Region Name : "Tally1"
c   Abaqus Part #      : 11
c   Abaqus Part Name   : "Tally1"
c   MCNP Pseudo-cell # : 11
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 21134.4 cm**3
c     # Cells          : 615
c     % of Total Cells : 0.69% 
c
c   Attila Region #    : 12
c   Attila Region Name : "Tally2"
c   Abaqus Part #      : 12
c   Abaqus Part Name   : "Tally2"
c   MCNP Pseudo-cell # : 12
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 35226 cm**3
c     # Cells          : 422
c     % of Total Cells : 0.48% 
c
c   Attila Region #    : 13
c   Attila Region Name : "Tally3"
c   Abaqus Part #      : 13
c   Abaqus Part Name   : "Tally3"
c   MCNP Pseudo-cell # : 13
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 49314.9 cm**3
c     # Cells          : 471
c     % of Total Cells : 0.53% 
c
c   Attila Region #    : 14
c   Attila Region Name : "Tally4"
c   Abaqus Part #      : 14
c   Abaqus Part Name   : "Tally4"
c   MCNP Pseudo-cell # : 14
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 200219 cm**3
c     # Cells          : 487
c     % of Total Cells : 0.55% 
c
c   Attila Region #    : 15
c   Attila Region Name : "Domain"
c   Abaqus Part #      : 15
c   Abaqus Part Name   : "Domain"
c   MCNP Pseudo-cell # : 15
c   Material           : "VOID_0"
c   Mesh Data
c     Meshed Volume    : 4.10785e+07 cm**3
c     # Cells          : 11806
c     % of Total Cells : 13.33% 
c
c ------------------------ End Input Information ---------------------------- 80
c
c ----------------------------- Cell Cards ---------------------------------- 80
1     0                  0                               u=1
2     1     -8           0                               u=1
3     1     -8           0                               u=1
4     1     -8           0                               u=1
5     0                  0                               u=1
6     0                  0                               u=1
7     0                  0                               u=1
8     0                  0                               u=1
9     0                  0                               u=1
10    0                  0                               u=1
11    0                  0                               u=1
12    0                  0                               u=1
13    0                  0                               u=1
14    0                  0                               u=1
15    0                  0                               u=1
16    0                  0                               u=1 $ background
17    0                  100 -101 102 -103 104 -105   fill=1 $ fill cell
18    0                  (-100:101:-102:103:-104:105)
c --------------------------- End Cell Cards -------------------------------- 80

c ---------------------------- Surface Cards -------------------------------- 80
c
100 px -155
101 px 155
102 py -155
103 py 155
104 pz -5
105 pz 705
c -------------------------- End Surface Cards ------------------------------ 80

c ----------------------------- Data Cards ---------------------------------- 80
c Embedded Geometry Specification
embed1 meshgeo=abaqus mgeoin=iter_port_plug_benchmark.abaq
       meeout=iter_port_plug_fwcadis.mcnp.eeout
       filetype=ascii
       background=16
       matcell= 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 &
       14 14 15 15
c
c Materials
c
c  Material 1: "Steel-el_1"
c  Constituents (weight %):
c fe-26000 (0.65115) ni-28000 (0.1225) mn55-25055 (0.018) mo-42000 (0.025)
c cu-29000 (0.003) cr-24000 (0.175) ti-22000 (0.00015) si-14000 (0.005)
m1  26000 -0.65115 28000 -0.1225 25055 -0.018 42000 -0.025 29000 -0.003 &
24000 -0.175 22000 -0.00015 14000 -0.005 
c
c Mode (Only n and/or p Currently Accepted)
mode n
c
c Cell Importances
imp:n 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
c
c Source Definition
c Written by Attila for MCNP 6.2.0 on:
c 15:54:31  23 Sep 2022 -04:00Z
c
sdef par=n
     pos=volumer
     erg=d1
     wgt= 1.0000000000E+00
c
c 1 source distribution(s) used to define the source.
c
si1 L
       1.4100000000E+01
sp1 D
       1.0000000000E+00
sb1 D
       1.0000000000E+00
c
c
c Histories (or Computer Time Cutoff)
nps 1
c ctme 1
c
c Imports WWP card(s) written by Attila:
read file=iter_port_plug_fwcadis.fwcadis.wwp
c
c
c Tallies or embee cards
c [UD]
c
c L'Ecuyer 63-bit random number generator (period=9.2E18)
rand gen=2
c
print -85 -86 -87 -98
c
c --------------------------- End Data Cards -------------------------------- 80
c End MCNP Input
