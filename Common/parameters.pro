//--------------------------------------------------------------------------------------------------------------
// Onelab Server Parameters/Variables:
//   These parameters enable easy interaction between GetDP/Gmsh (ONELAB)
//   and python scripts 
//--------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------
// Geometry parameters
//--------------------------------------------------------------------------------------------------------------
inputname_R_c     = "Input/0Geometry Parameters/2Conductor Radius (R_c)/(m)";
inputname_rs      = "Input/0Geometry Parameters/3Sleeve Radius (rs)/(m)";
inputname_distx   = "Input/0Geometry Parameters/4Distance Between Conductors in X/(m)";
inputname_disty   = "Input/0Geometry Parameters/5Distance Between Conductors in Y/(m)";
inputname_rows    = "Input/0Geometry Parameters/7Coil Bundle Array/0Rows";
inputname_columns = "Input/0Geometry Parameters/7Coil Bundle Array/1Columns"; 
inputname_nturns  = "Input/0Geometry Parameters/7Coil Bundle Array/2Number of Turns";

//-------------------------------------------------------------------------------------------------------------- 

//--------------------------------------------------------------------------------------------------------------
// Model Options:
//--------------------------------------------------------------------------------------------------------------
inputname_model       = "Input/0Model Configuration/0Type";
//-------------------------------------------------------------------------------------------------------------- 

//--------------------------------------------------------------------------------------------------------------
// Particular test cases: 2D/3D, Outline of the conductor and curve conductor
//--------------------------------------------------------------------------------------------------------------
inputname_condOutline = "Input/1Geometry Parameters/2Feature Checkboxes/Show Conductor Outline";
//--------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------
// Material Properties and Frequency parameters
//--------------------------------------------------------------------------------------------------------------
inputname_freq           = "Input/2Model Excitation/3Frequency (Hz)/ ";
inputname_Vp             = "Input/2Model Excitation/1Voltage (V)/Wire %g";
inputname_Ip             = "Input/2Model Excitation/2Current (A)/Wire %g";
inputname_Vg             = "Input/2Model Excitation/1Voltage Ground (V)/Wire %g";
inputname_excitation     = "Input/2Model Excitation/0Excitation Type";
inputname_connectiontype = "Input/2Model Excitation/4Feature Checkboxes/Connect in Series";

inputname_sig     = "Input/1Conductor Material/ Electrical Conductivity (S\m)/Wire %g";
inputname_mur     = "Input/1Conductor Material/ Relative Permeability/ ";
inputname_epr     = "Input/1Conductor Material/ Relative Permittivity/ ";

//--------------------------------------------------------------------------------------------------------------
// Output variables
//--------------------------------------------------------------------------------------------------------------
outputname_resistance_source  = "Output/0Impedance/0Resistance at Source (2D)/0Full";
outputname_inductance_source  = "Output/0Impedance/1Inductance at Source (2D)/0Full";

outputname_resistance_wire  = "Output/0Impedance/0Resistance (2D)/0Full";
outputname_inductance_wire  = "Output/0Impedance/1Inductance (2D)/0Full";

outputname_resistance_wire_thin  = "Output/0Impedance/0Resistance (2D)/1Corr";
outputname_inductance_wire_thin  = "Output/0Impedance/1Inductance (2D)/1Corr";

// 
//outputname_resistance_prox  = "}Resistance";
//outputname_rePower_prox = "}Real Power";
outputname_resistance_prox  = "Output/0Impedance/0Resistance (2D)/1Corr";
outputname_rePower_prox = "Output/0Impedance/1Inductance (2D)/1Corr";

 
// Physical Entities Tags

//2D
Full_DomainCC    = 1000000;
Full_DomainC     = 2000000;
Full_skinDomainC = 3000000;
LR_DomainC       = 4000000;
Struct_Sleeve    = 5000000;
Inf              = 6000000;
noFlux           = 7000000;


// Initial Parameters
columns0 = 3;
rows0 = 1;

disty0 = 8e-3;//2.05e-3;//8e-3; //12e-3
distx0 = disty0;
R_c0 =  1e-3;
rs0 = 1e-3;

//--------------------------------------------------------------------------------------------------------------
// Initial Material Properties
//--------------------------------------------------------------------------------------------------------------
  eps0 = 8.85e-12;          // Free Space Permittivity 
  mu0 = 4*Pi*1e-7;          // Free Space Permeability
  
  sig0  = 5.96e7;           // Conductivity of Copper
  mur0  = 1;                // Relative Permeability
  epr0  = 1;                // Relative Permittivity
  
DefineConstant[   sig    = {sig0  ,                    Name StrCat[inputname_sig]                        , Visible 0}, 
                  mur    = {mur0  ,                    Name StrCat[inputname_mur]                        , Visible 0},
		  epr    = {epr0  ,                    Name StrCat[inputname_epr]                        , Visible 0}
		];
  
//--------------------------------------------------------------------------------------------------------------
//Initial Global Sources (Peak Values)  and Frequency
//--------------------------------------------------------------------------------------------------------------

  //Sources
  V0    = 1;              // Peak Voltage
  I0    = 1;              // Peak Current
  Freq0 = 100000;              // Initial Frequency of the problem

//--------------------------------------------------------------------------------------------------------------

air_distance = 40e-3;

// Static Menus
isCurved = 0;
isRound = 1;
isBoundaryLayer = 1;
withIED = 0 ;
withGlobal = 1;
AirAxi = 1;

 //--------------------------------------------------------------------------------------------------------------
// Geometrical Menu
//--------------------------------------------------------------------------------------------------------------
   DefineConstant[ 
		  R_c     = {R_c0,     Min R_c0/10,   Max 10*R_c0    , Step 1, Name StrCat[inputname_R_c]                            , Visible 0},
		  distx   = {distx0,   Min distx0,    Max 1000*distx0, Step 1, Name StrCat[inputname_distx]                          , Visible 1},
		  disty   = {disty0,   Min disty0/10, Max 1000*disty0, Step 1, Name StrCat[inputname_disty]                          , Visible 1},
		  rows    = {rows0,    Min 1,         Max 15,          Step 1, Name StrCat[inputname_rows]   , Highlight "AliceBlue" , Visible 1},
		  columns = {columns0, Min 1,         Max 10,          Step 1, Name StrCat[inputname_columns], Highlight "AliceBlue" , Visible 1},
		  Freq   = {Freq0 , Min 1e-3, Max 1e9, Name StrCat[inputname_freq], Highlight "AliceBlue", Visible 1}
		   ];
 
// Calculate and show number of turns		   
DefineConstant[NumWires    = {rows*columns, Name StrCat[inputname_nturns], Highlight "LightGreen" , Visible 1, ReadOnly}];
 
 //--------------------------------------------------------------------------------------------------------------
// Modeling Approach: 
//                    Full Model (Original Volume Model) = 0
//                    Thin approximation without correction = 1 
//                    Thin Approximation with correction accounting with skin and proximity effects = 2
//--------------------------------------------------------------------------------------------------------------

DefineConstant[
	       isThin = {2, Choices{
		         0="Full_Model",
		         1="LR_NoCorrection",
			 2="LR_SemiAnalytical"},
    Name StrCat[inputname_model] , Highlight "Blue", Visible 1}  	 
  ];

  

If (isThin)
    DefineConstant[
        withRing = {1, Choices {0,1},
            Name StrCat[inputname_condOutline], Visible  1 }
       
 ];
 
     DefineConstant[ rs = {rs0, Min rs0, Max 3*rs0 , Step 1, 
	Name StrCat[inputname_rs], Visible withRing}];
Else
    DefineConstant[
        withRing = {1, Choices {0,1},
            Name StrCat[inputname_condOutline], Visible  0 }
 ];
 
     DefineConstant[ rs = {rs0, Min rs0, Max 3*rs0 , Step 1, 
	Name StrCat[inputname_rs], Visible 0}];
EndIf


//--------------------------------------------------------------------------------------------------------------
// Excitation Source Type
//--------------------------------------------------------------------------------------------------------------
		
  DefineConstant[
	       sourceType = {0, Choices{
		         0="Current",
		         1="Voltage"},
    Name StrCat[inputname_excitation] , Highlight "Blue", Visible 1}  	 
  ];	
  

//--------------------------------------------------------------------------------------------------------------
// Connect all conductors in series as turns (single multiturn coil) to a circuit source
//--------------------------------------------------------------------------------------------------------------
 DefineConstant[
        isSeries = {1, Choices {0,1},
            Name StrCat[inputname_connectiontype], Visible  1 }
    ]; 

    
    
    
//--------------------------------------------------------------------------------------------------------------
// If conductors are in parallel identify which conductors are active
//--------------------------------------------------------------------------------------------------------------
active_wire = 1;

//flush every excitation menu (maximum 200 for now)
For i In {NumWires+1:200}
  DefineConstant[ Ip~{i} = {0  , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, i], Visible 0}];
  DefineConstant[ Vp~{i} = {0  , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, i], Visible 0}];
EndFor

If(!isSeries)
    If(sourceType == 0)
        For i In {1:NumWires}
            If(i != active_wire)
	        param = 0;
	    Else
	        param = 1;
	    EndIf
            DefineConstant[ Ip~{i} = {I0*param   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, i], Visible 1}];
	    DefineConstant[ Vp~{i} = {V0*param   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, i], Visible 0}];
            ActiveSource~{i} =  Ip~{i};
        EndFor
    ElseIf(sourceType == 1)
        For i In {1:NumWires}
            If(i != active_wire)
	        param = 0;
	    Else
	        param = 1;
	    EndIf
            DefineConstant[ Vp~{i} = {V0*param   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, i], Visible 1}];
            DefineConstant[ Ip~{i} = {I0*param   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, i], Visible 0}];
	    ActiveSource~{i} =  Vp~{i};
        EndFor
    EndIf

Else
      
For i In {2:200}
  DefineConstant[ Ip~{i} = {0  , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, i], Visible 0}];
  DefineConstant[ Vp~{i} = {0  , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, i], Visible 0}];
EndFor

    //Note that this is temporary. I'm assuming that it's a single coil
    If(sourceType == 0)
          DefineConstant[ Ip~{1} = {I0   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, 1], Visible 1}];
	  DefineConstant[ Vp~{1} = {V0    , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, 1], Visible 0}];
          ActiveSource~{1} =  Ip~{1};
     ElseIf(sourceType == 1) 
          DefineConstant[ Vp~{1} = {V0   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vp, 1], Visible 1}];
          DefineConstant[ Ip~{1} = {I0   , Min 1e-3, Max 1e9, Name Sprintf[inputname_Ip, 1], Visible 0}];
	  ActiveSource~{1} =  Vp~{1};
    EndIf
    
EndIf

For i In {1:NumWires}
    DefineConstant[ Vg~{i} = {0    , Min 1e-3, Max 1e9, Name Sprintf[inputname_Vg, i], Visible 0}];
EndFor


    
//--------------------------------------------------------------------------------------------------------------
// Postprocessing 
//--------------------------------------------------------------------------------------------------------------


// Show color maps
DefineConstant[
    showMaps = {1, Choices {0,1},
        Name "PostProcessing/Feature Checkboxes/Show Maps", Visible  1 }
];


//--------------------------------------------------------------------------------------------------------------
// Setup of centerlines and distances between conductors or conductor turns
//--------------------------------------------------------------------------------------------------------------
If(rows == 0)
    dist_row = 0;
Else
    dist_row = -(rows-1)*(disty/2);
EndIf

dist_row_initialize = dist_row;

If(AirAxi)
    dist_col = 0.005;
Else
    dist_col = -(columns-1)*(distx/2);
EndIf

// Initial separation row/column are the bottom and top of line plot
bottom_y = dist_row - R_c;
top_y = -dist_row + R_c;
left_x = dist_col;
right_x = dist_col;

// Creation of conductor center-line coordinates
For j In {1:columns} 
    For i In {1:rows} 
                
	cx~{i}~{j} = 0 + dist_col;
	cy~{i}~{j} = 0 + dist_row;
	
	// create distance between each row 
	If(j%2 == 0)
            dist_row -= disty;
	Else
	    dist_row += disty;
        EndIf
	
    EndFor
    // re-initialize row_distance
    If(j%2 == 0)
        dist_row = dist_row_initialize;
    Else
        dist_row = -dist_row_initialize;
    EndIf
    
    // create distance between each column
    dist_col += distx;   
    
EndFor

//-----

/*
Calculation of skin depth and 
meshing parameter based on frequency variation
*/
omega = 2*Pi*Freq;
skin_depth =  1/((Sqrt[2]/2) * Sqrt[omega*sig*(mu0*mur)]);

If(skin_depth >= R_c)
    mesh_parameter = R_c/3;
Else
    mesh_parameter = skin_depth/3; // 3 elements per skin depth
EndIf

// Show color maps
DefineConstant[
    freqRefinement = {1, Choices {0,1},
        Name "PreProcessing/Feature Checkboxes/Frequency Mesh Refinement", Visible  1 }
];
