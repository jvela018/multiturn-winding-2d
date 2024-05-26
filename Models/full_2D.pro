//--------------------------------------------------------------------------------------------------------------
// Full 2D Model/Mesh
// EM Formulations: Magnetodynamics
// Simple example on series connection
//
//--------------------------------------------------------------------------------------------------------------


// Create Folders for Results
 DefineConstant[FileId = ""];
DefineConstant[ResultsDir  = StrCat["getdp_results_full_2D/",FileId]];

Group{
   // Create independent regions for each conductor
    For i In {0:NumWires-1}
        //Full Model Domains -----------------------------------------------------------------------
        // Create entities for each individual coil
        Full_DomainC~{i+1}  += Region[{(Full_DomainC+i)}];
        Full_Anode~{i+1}    += Region[{(Full_DomainC+i)}];

        // Group each entity into classic conducting and non-conducting domain DomainC and DomainCC
        DomainC  += Region[{Full_DomainC~{i+1}}];
        Anodes   += Region[{Full_Anode~{i+1}}];
	Electrodes_with_GQ += Region[{Full_Anode~{i+1}}];
    EndFor

  
  If(withIED)
      IED = Region[{IED}];
      DomainCC += Region[{IED,Full_DomainCC}];  
  Else
      DomainCC += Region[{Full_DomainCC}];  
  EndIf
  
  
  Domain   = Region[{DomainC, DomainCC}];
  Infinity = Region[{Inf}];
  noFlux = Region[{noFlux}];


  

  // Circuit network groups

  // Sources
  CurrentSource1 = Region[{3001}];
  VoltageSource1 = Region[{3002}];

  If(sourceType == 0)
    Sources_Cir = Region[{CurrentSource1}];
  ElseIf(sourceType == 1)
    Sources_Cir = Region[{VoltageSource1}];
  EndIf

  // Complete circuit domain containing all circuit elements
  Domain_Cir = Region[{Sources_Cir}];
}


//--------------------------------------------------------------------------------------------------------------
// Functions
//--------------------------------------------------------------------------------------------------------------

Function{
//--------------------------------------------------------------------------------------------------------------
// Material Properties
//--------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------
// Initial Material Properties
//--------------------------------------------------------------------------------------------------------------
  eps0 = 8.85e-12;             // Free Space Permittivity
  mu0 = 4*Pi*1e-7;             // Free Space Permeability

  sig  = 5.96e7;               // Conductivity of Copper
  mur  = 1;                    // Relative Permeability
//--------------------------------------------------------------------------------------------------------------

  sigma[DomainC]   = sig;      //Electrical Conductivity
  sigma[DomainCC]  = 0.;       //Electrical Conductivity in Non-Conducting Domain

  mu[DomainC]      = mur*mu0;  //Permeability of Free Space
  mu[DomainCC]     = mur*mu0;  //Permeability of Conducting Domain
  nu[]             = 1/mu[];   //Reluctivity

//--------------------------------------------------------------------------------------------------------------
// Excitation Parameters
//--------------------------------------------------------------------------------------------------------------
 i[] = Complex[0,1];
 //omega = 2*Pi*Freq;
 //skin_depth =  1/((Sqrt[2]/2) * Sqrt[omega*sig*(mu0*mur)]);
  
  DefineConstant[ skind   = {skin_depth, Name "Input/Skin Effect/0Skin Depth (m) " , Highlight "LightGreen", Visible 1, ReadOnly}
                  condrad = {R_c  , Name "Input/Skin Effect/1Conductor Radius (m) ", Highlight "LightGreen", Visible 1, ReadOnly}
		  skin_rad_ratio = {R_c/skin_depth  , Name "Input/Skin Effect/2Ratio", Highlight "LightGreen", Visible 1, ReadOnly}
		  mesh_para = {mesh_parameter  , Name "Input/Skin Effect/3 Mesh Characteristic Length", Highlight "LightGreen", Visible 1, ReadOnly}
  ];
  
  For i In {1:NumWires}
      If(isSeries)
          Vin~{i}[] = I0;
          Iin~{i}[] = V0;
      Else
	  Vin~{i}[] = Vp~{i};
          Iin~{i}[] = Ip~{i};  
      EndIf
  EndFor
  
  
  // Function to extract degrees of freedom and print them on separate file (only used in resolution)
   myRegionNumberList() = GetRegions[Domain];
   nElements[] = GetNumElements[]{myRegionNumberList()};
}


//--------------------------------------------------------------------------------------------------------------
// Jacobian: Note that the Jacobian is set up for Full Model Integrals
//--------------------------------------------------------------------------------------------------------------

Jacobian {
  { Name Vol ;
    Case {      

      If(isCurved && withIED)                                             // isCurved will model a sphere in 2D
          { Region All ; Jacobian VolAxiSqu ; } 
	  { Region IED ; Jacobian VolAxiSphShell {Rint,Rext} ; }        
      ElseIf(isCurved && !withIED)
	  { Region All ; Jacobian VolAxiSqu ; }
	ElseIf(!isCurved && !withIED)                                      // !isCurved models a cylinder  in 2D
	  { Region All ; Jacobian Vol ; }	  
      ElseIf(!isCurved && withIED)
	  { Region All ; Jacobian Vol ; }
          { Region IED ; Jacobian VolSphShell {Rint,Rext} ; }        
      EndIf
    }
  }
}


//--------------------------------------------------------------------------------------------------------------
// Integration
//--------------------------------------------------------------------------------------------------------------
Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
			{ GeoElement Point       ; NumberOfPoints  1; }
		        { GeoElement Line        ; NumberOfPoints  3; }
                        { GeoElement Triangle    ; NumberOfPoints  3; }
			{ GeoElement Quadrangle  ; NumberOfPoints  4 ; }
			{ GeoElement Prism       ; NumberOfPoints  21 ; }
			{ GeoElement Tetrahedron ; NumberOfPoints  4 ; }
                        { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
                        { GeoElement Pyramid     ; NumberOfPoints  8 ; }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------------------
// Constraints (Dirichlet Boundary Conditions)
//--------------------------------------------------------------------------------------------------------------

Constraint{


  { Name MVP_2D ;
    Case {
      { Region Infinity ; Type Assign ; Value 0. ; }
      { Region noFlux ; Type Assign ; Value 0. ; }
      
    }
  }

  { Name Current_2D ;
    Case {
      If(!isSeries && sourceType == 0)
          For i In {1:NumWires}
	      { Region Full_Anode~{i};    Value Iin~{i}[]; }
          EndFor
      EndIf
    }
  }

  { Name Voltage_2D ;
    Case {
      If(!isSeries && sourceType == 1)
          For i In {1:NumWires}
	      { Region Full_Anode~{i};    Value Vin~{i}[]; }
          EndFor
      EndIf
    }
  }


  // Constraints related to circuit connections

  { Name Current_Cir ;
    Case {
      If(isSeries && sourceType == 0)
	{ Region CurrentSource1; Value -I0; }
      EndIf
    }
  }
  { Name Voltage_Cir ;
    Case {
      If(isSeries && sourceType == 1)
	{ Region VoltageSource1; Value V0; }
      EndIf
    }
  }

  
  If(isSeries)
  { Name ElectricalCircuit ; Type Network ;
    Case Circuit1 {
      
      If(sourceType == 1)
        { Region VoltageSource1;   Branch{1, 2}; }
      ElseIf(sourceType == 0)
        { Region CurrentSource1;   Branch{2, 1}; }
      EndIf
      
     For i In {1:NumWires}
          node_1 = i+1;
          node_2 = i+2;
	  If (i == NumWires)
	     { Region Full_Anode~{i};     Branch{node_1, 1}; }
	  Else
	    { Region Full_Anode~{i};     Branch{node_1, node_2}; } 	
	  EndIf
      EndFor
    
    }
  }
 EndIf

}

//--------------------------------------------------------------------------------------------------------------
// 2D Function Spaces
//--------------------------------------------------------------------------------------------------------------


FunctionSpace {

    // Magnetic Vector Potential
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Region[{Domain}] ; Entity NodesOf [ All ] ; }
   }
    Constraint {
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }

  // Gradient of Electric scalar potential (2D)
  { Name Hregion_u_2D ; Type Form1P ;
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC ; Entity DomainC ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType GroupsOfNodesOf ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ; EntityType GroupsOfNodesOf ; NameOfConstraint Current_2D ; }
    }
  }

  // Circuit related function space
   { Name Hregion_Cir; Type Scalar;
    BasisFunction {
      { Name sn; NameOfCoef ir; Function BF_Region;
        Support Domain_Cir; Entity Domain_Cir; }
    }
    GlobalQuantity {
      { Name Iz; Type AliasOf; NameOfCoef ir; }
      { Name Uz; Type AssociatedWith; NameOfCoef ir; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

}


//--------------------------------------------------------------------------------------------------------------
// Electromagnetic (EM) Formulations: Magnetodynamics
//--------------------------------------------------------------------------------------------------------------

Formulation {

  { Name EM_Formulations; Type FemEquation;
    Quantity {

      { Name a; Type Local; NameOfSpace Hcurl_a_2D; }
      { Name v; Type Local; NameOfSpace Hregion_u_2D; }
      { Name U; Type Global; NameOfSpace Hregion_u_2D [U]; }
      { Name I; Type Global; NameOfSpace Hregion_u_2D [I]; }

      //Circuit related info
      If(isSeries)
          { Name Uz; Type Global; NameOfSpace Hregion_Cir [Uz]; }
          { Name Iz; Type Global; NameOfSpace Hregion_Cir [Iz]; }
      EndIf
    }

    Equation {

       Galerkin { [ nu[] * Dof{d a} , {d a} ];
        In Domain; Jacobian Vol; Integration I1; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ];
        In DomainC; Jacobian Vol; Integration I1; }
      Galerkin { [ sigma[] * Dof{v} , {a} ];
        In DomainC; Jacobian Vol; Integration I1; }
      Galerkin { DtDof[sigma[] * Dof{a} , {v} ];
        In DomainC; Jacobian Vol; Integration I1; }


	 Galerkin { [ sigma[] * Dof{v} , {v} ];
        In DomainC; Jacobian Vol; Integration I1; }

	GlobalTerm { [ Dof{I}, {U} ]; In Electrodes_with_GQ; }

        If(isSeries)
          //+++ No meaning      GlobalTerm { [ Dof{Iz} , {Uz} ]; In currentSource; }
          //+++ Needed now for Dof{Iz} to be defined (no effect because multiplied by 0)
          GlobalTerm { [ 0 * Dof{Iz} , {Iz} ]; In Domain_Cir; }

          // Here is the connection between the circuit and the finite element simulation
          GlobalEquation { Type Network; NameOfConstraint ElectricalCircuit;
            { Node {I};  Loop {U};  Equation {I};  In Electrodes_with_GQ; }
            { Node {Iz}; Loop {Uz}; Equation {Uz}; In Domain_Cir; }
          }
        EndIf


    }
  }
}
//--------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------
// Resolutions
//--------------------------------------------------------------------------------------------------------------

  Resolution {

 { Name Resolutions;
    System {
      { Name S; NameOfFormulation EM_Formulations;
	Type ComplexValue; Frequency Freq;}
    }
    Operation {
      CreateDir[StrCat['../',ResultsDir]];
      SetGlobalSolverOptions["-petsc_prealloc 2000"];
      
      Generate[S]; Solve[S]; 
      
      Print[{Freq, $KSPSystemSize, nElements[]}, Format "%g %g %g", File StrCat[ResultsDir, "2D_Dof_nElem_Full.dat"] ]; 
      
      SaveSolution[S]; 
      
      // Automatic Postprocessing
      For i In {1:NumWires}
	  PostOperation[ConductorArea~{i}];
	  PostOperation[Flux~{i}];
      EndFor
      
      // Command Window Prints
      For i In {1:NumWires}
      Print[ {i,$conductor_area~{i},Sqrt[$conductor_area~{i}/Pi], $flux~{i}, disty  }, 
                 Format "Wire %g: conductor_area = %7.5e,conductor_radius = %7.5e, flux = %7.5e , distance = %7.5e"] ;
      EndFor
      
      PostOperation[PostOperations];
    }
  }
}


//--------------------------------------------------------------------------------------------------------------
// Postprocessing
//--------------------------------------------------------------------------------------------------------------

PostProcessing {

  { Name PostProcessings; NameOfFormulation EM_Formulations;
    PostQuantity {


//--------------------------------------------------------------------------------------------------------------
// Local Quantities:
//              v: Electric Scalar Potential (V)
//              a: Magnetic Vector Potential (W/m)
//              b: Magnetic Flux Density  (T)
//              jz: z component of current Density (A/m^2)
//--------------------------------------------------------------------------------------------------------------
    { Name j;  Value { Term { [ sigma[]*(Dt[{a}] + {v}) ]; In DomainC ; Jacobian Vol; } } }
    { Name v;   Value { Term { [ {v} ];                     In DomainC; Jacobian Vol; } } }
    { Name b;   Value { Term { [ {d a} ];                   In Domain ; Jacobian Vol; } } }
    { Name az;  Value { Term { [ CompZ[{a}] ];              In Domain ; Jacobian Vol; } } }


//--------------------------------------------------------------------------------------------------------------
// Global Quantities:
//              reU: Real component of Voltage (V)
//              reI: Real component of Current (A)
//              imU: Imaginary component of Voltage (V)
//              imI: Imaginary component of Current (A)
//--------------------------------------------------------------------------------------------------------------

     { Name reU;  Value { Term { [ Re[{U}] ];  In Electrodes_with_GQ; } } }
     { Name reI;  Value { Term { [ Re[{I}] ];  In Electrodes_with_GQ; } } }
     { Name imU;  Value { Term { [ Im[{U}] ];  In Electrodes_with_GQ; } } }
     { Name imI;  Value { Term { [ Im[{I}] ];  In Electrodes_with_GQ; } } }
     
     // Retreives the resistance, self and mututual inductance in the case of parallel wires (inductance matrix)
     If(!isSeries)
         For j In {1:NumWires}	  
	     If(ActiveSource~{j} !=0 )		
	         For i In {1:NumWires}  
		     If(i==j)
		       
	                 { Name R~{i}~{j}; Value { Term { [-Re[{U}/{I}]];        In Full_Anode~{j}; } } }
		         { Name L~{i}~{j}; Value { Term { [-Im[{U}/{I}]/omega];  In Full_Anode~{j}; } } }
		    
	             Else
		      
		         { Name R~{i}~{j}; Value { Term { [Re[i[]*omega*$flux~{i}/{I}]];        In Full_Anode~{j}; } } }
                         { Name L~{i}~{j}; Value { Term { [ Im[i[]*omega*$flux~{i}/{I}]/omega]; In Full_Anode~{j}; } } }                  
		      
		     EndIf
                 EndFor	    
	     EndIf
	 EndFor	  
     EndIf  
  
     
     If(isSeries)
         
         { Name Rz;  Value { Term { [ Re[{Uz}/{Iz}] ];        In Sources_Cir; } } }
         { Name Lz;  Value { Term { [ Im[{Uz}/{Iz}]/omega ];  In Sources_Cir; } } }
       
	 
         { Name reUz;  Value { Term { [ Re[{Uz}] ];  In Sources_Cir; } } }
         { Name reIz;  Value { Term { [ Re[{Iz}] ];  In Sources_Cir; } } }
         { Name imUz;  Value { Term { [ Im[{Uz}] ];  In Sources_Cir; } } }
         { Name imIz;  Value { Term { [ Im[{Iz}] ];  In Sources_Cir; } } }
	 
     EndIf
     
     // Separate PostOperations
     For i In {1:NumWires}   
        { Name getArea~{i}; Value{ Integral{ [ 1 ] ; In Full_Anode~{i}; Integration I1 ; Jacobian Vol; } } }
        { Name getFlux~{i}; Value { Integral {[ CompZ[{a}]/$conductor_area~{i} ] ; In Full_Anode~{i}; Jacobian Vol ; Integration I1 ; } } }
     EndFor

     }
   }
 }

//--------------------------------------------------------------------------------------------------------------
// PostOperation
//--------------------------------------------------------------------------------------------------------------


PostOperation PostOperations UsingPost PostProcessings {
      
  If(showMaps)
       Print [j    , OnElementsOf DomainC, File  StrCat[ResultsDir,"j_MagDyn.pos"] ];
       Print [v    , OnElementsOf DomainC, File StrCat[ResultsDir,"v_MagDyn.pos"]  ];
       Print [b    , OnElementsOf Domain,  File StrCat[ResultsDir,"b_MagDyn.pos"]];
       
       
       Print [az   , OnElementsOf Domain,  File  StrCat[ResultsDir,"az_MagDyn.pos"] ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;",
		 "View[l].IntervalsType = 3;",
		 "View[l].NbIso = 30;",
		 "View[l].NormalRaise = 0;",
	         "View[l].RangeType = 3;"],
       File  StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly] ;
     
     If(isSeries)
       Print [ az, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
	  Format Gmsh, File StrCat[ResultsDir,"Cut_az.pos"] ];
	
	 Print [ az, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
	  Format SimpleTable, File StrCat[ResultsDir,"Cut_az.dat"] ];
	
	// Some bs test
	 Print [ az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	  Format Gmsh, File StrCat[ResultsDir,"Cut_az.pos"] ];
	
     Else
	  Print [ az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	  Format Gmsh, File StrCat[ResultsDir,"Cut_az.pos"] ];
     EndIf
	
		Echo[ StrCat["l=PostProcessing.NbViews-1;",
	       "View[l].Axes = 3;",
	       "View[l].LineWidth = 3;",
	       "View[l].Type = 2;"],
	       File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly];
	     
	 	 Print [  az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	  Format SimpleTable, File StrCat[ResultsDir,"Cut_az.dat"] ];
     
  EndIf

       
	
	// Total voltage and current measured from source
	If(isSeries)
	    Print [Rz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Rz.dat"], 
	      SendToServer StrCat[outputname_resistance_source], Color "AliceBlue"];
	    
	    Print [Lz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Lz.dat"], 
	      SendToServer StrCat[outputname_inductance_source], Color "AliceBlue"];
	  Else

	      // Retreive resistance, self and mutual inductance
	      For j In {1:NumWires}
                  If(ActiveSource~{j} !=0 )
	              For i In {1:NumWires}
		       
	                  ext_ij = StrCat[Sprintf["%g",i],Sprintf["%g",j]]; // Zji 
	                  If(i==j)
	                      Print[ R~{i}~{j}, OnRegion  Full_Anode~{j}, File > StrCat[ResultsDir,"Rper_full.dat"], Format Table,
		                  SendToServer StrCat[outputname_resistance_wire, ext_ij] ];		   

	                  EndIf
	       
	                      Print[ L~{i}~{j}, OnRegion  Full_Anode~{j}, File > StrCat[ResultsDir,"Lper_full.dat"], Format Table,
		                  SendToServer StrCat[outputname_inductance_wire, ext_ij] ];
		       
	              EndFor
	          EndIf
	      EndFor
	     
	 EndIf
	
 	
     
}


For i In {1:NumWires} 

 PostOperation ConductorArea~{i} UsingPost PostProcessings {
     Print[ getArea~{i}[Full_Anode~{i}], OnGlobal, 
        Format Table, File "", StoreInVariable $conductor_area~{i}];
}
  
 PostOperation Flux~{i} UsingPost PostProcessings {
  Print[ getFlux~{i}[Full_Anode~{i}], OnGlobal, 
    Format Table, File "", StoreInVariable $flux~{i}];
 }
EndFor


DefineConstant[
  R_ = {"Resolutions"    , Name "GetDP/1ResolutionChoices"      , Visible 1},
  C_ = {"-solve -v2"     , Name "GetDP/9ComputeCommand"         , Visible 1},
  P_ = {""               , Name "GetDP/2PostOperationChoices"   , Visible 1}
];


