//--------------------------------------------------------------------------------------------------------------
// Full 2D Model/Mesh
// EM Formulations: Magnetodynamics
// Simple example on series connection
//
//--------------------------------------------------------------------------------------------------------------


// Create Folders for Results
DefineConstant[FileId = ""];
DefineConstant[ResultsDir  = StrCat["getdp_results_thin_mat_2D/",FileId]];


Group{

   // Create independent regions for each conductor
    For i In {0:NumWires-1}
        //Full Model Domains -----------------------------------------------------------------------
        // Create entities for each individual coil
	
	If(withRing)
            Full_DomainC~{i+1}  += Region[{(Full_DomainC+i)}];
	Else
	    Full_DomainC~{i+1}  += Region[{}];  
	EndIf
	LR_DomainC~{i+1}  += Region[{(LR_DomainC+i)}];
        LR_Anode~{i+1}    += Region[{(LR_DomainC+i)}];

        // Group each entity into classic conducting and non-conducting domain DomainC and DomainCC
        DomainC  += Region[{LR_DomainC~{i+1}}];
	LR_DomainC += Region[{LR_DomainC~{i+1}}];
        Anodes   += Region[{LR_Anode~{i+1}}];
	Electrodes_with_GQ += Region[{LR_Anode~{i+1}}];
	voidRegion  += Region[{Full_DomainC~{i+1}}];
    EndFor

  
  If(withIED)
      IED = Region[{IED}];
      DomainCC += Region[{IED,Full_DomainCC}];  
  Else
      DomainCC += Region[{Full_DomainCC}];  
  EndIf
  
  
  Domain   = Region[{DomainC, DomainCC, voidRegion}];
  noFlux = Region[{noFlux}];
  Infinity = Region[{Inf}];

  
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
  
  // support
  Dom_Hcurl_a = Region[ Domain ];
  Dom_Hthin_a = ElementsOf[ Domain, OnOneSideOf LR_DomainC ]; // Sleeve
  Dom_Hregion_i = Region[ LR_DomainC ];
  
   //Sleeve Domain Per Turn
  For i In {1:NumWires}
      Dom_Sleeve~{i} = ElementsOf[ voidRegion, OnOneSideOf LR_DomainC~{i} ];
  EndFor
  
  //All Sleeves
  Dom_Sleeve = ElementsOf[ DomainCC, OnOneSideOf LR_DomainC ];
        
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

//sigma[DomainC]   = sig;      //Electrical Conductivity
  sigma[DomainCC]  = 0.;       //Electrical Conductivity in Non-Conducting Domain

  //mu[DomainC]      = mur*mu0;  //Permeability of Free Space
  mu[Region[{DomainCC,voidRegion}]]     = mur*mu0;  //Permeability of Conducting Domain
  nu[]             = 1/mu[];   //Reluctivity

//--------------------------------------------------------------------------------------------------------------
// Excitation Parameters
//--------------------------------------------------------------------------------------------------------------

// Complex unit i
  i[] = Complex[0,1];
 //Area of the conductors
  A_s = Pi*rs^2;
  A_c = Pi*R_c^2;
  
  // Resistance and Internal Inductance due to skin effect
  R_dc = 1/(sig*A_c);
  //omega = 2*Pi*Freq;
  //skin_depth =  1/((Sqrt[2]/2) * Sqrt[omega*sig*(mu0*mur)]);
  
  k[] = (1-i[])/skin_depth;
  k_conj[] = Conj[k[]];

  //Bessel Function definition
  J_0[] = JnComplex[0,$1];
  J_1[] = JnComplex[1,$1];
  J_2[] = JnComplex[2,$1];
  //-------------------------
  
  //Effective Materials
  j0[] = I0/(A_c); 
  js[] = (k[]*I0)/(2*Pi*R_c)*(J_0[k[]*R_c]/J_1[k[]*R_c]);
  sigma[DomainC]   = sig*(j0[]/js[]);      //Electrical Conductivity
  
  nu_effective[] = ((1/mu0)*(1 + mu0*((2*i[]*omega*sig*R_c^3*J_2[k[]*R_c])/(R_c*k[]^2*(J_0[k[]*R_c] - J_2[k[]*R_c])+2*k[]*J_1[k[]*R_c]))*A_s))/A_c;  // Magnetic Reluctivity 
  //mur_effective[] = mur*((J_1[k[]*R_c])/(J_0[k[]*R_c]-J_1[k[]*R_c]));
  //nu_effective[] = 1/(mu0*mur_effective[])/A_s;// Magnetic Reluctivity 

  
  DefineConstant[ skind   = {skin_depth, Name "Input/Skin Effect/0Skin Depth (m) " , Highlight "LightGreen", Visible 1, ReadOnly}
                  condrad = {R_c  , Name "Input/Skin Effect/1Conductor Radius (m) ", Highlight "LightGreen", Visible 1, ReadOnly}
		  skin_rad_ratio = {R_c/skin_depth  , Name "Input/Skin Effect/2Ratio", Highlight "LightGreen", Visible 1, ReadOnly}
		  mesh_para = {mesh_parameter  , Name "Input/Skin Effect/3 Mesh Characteristic Length", Highlight "LightGreen", Visible 0, ReadOnly}
  ];
  
  Vin[] = V0;
  Iin[] = I0;
  


  CoefJac[] = A_c;
  
  For i In {1:NumWires}
      If(isSeries)
          Vin~{i}[] = I0;
          Iin~{i}[] = V0;
      Else
	  Vin~{i}[] = Vp~{i};
          Iin~{i}[] = Ip~{i};  
      EndIf
  EndFor
  
  //Prox stuff from Joule Loss
  p_denom[] = (R_c*k_conj[]*(J_0[k_conj[]*R_c] - J_2[k_conj[]*R_c]) + 2*J_1[k_conj[]*R_c])*(R_c*k[]*(J_0[k[]*R_c] - J_2[k[]*R_c]) + 2*J_1[k[]*R_c]) * (k[]^2 - k_conj[]^2);
  p_nume[] = (16*omega^2*sig*R_c^2)*(Pi*R_c*(k_conj[]*J_1[k[]*R_c]*J_0[k_conj[]*R_c] - k[]*J_0[k[]*R_c]*J_1[k_conj[]*R_c]));
  p_prox[] = p_nume[]/p_denom[]; 
  
  
  
  
// cylindrical r coordinate for correction of magnetic vector potential in post-processing 
n = 1;
For i In {1:rows}
    For j In {1:columns}
        // converting x and y coordinates to r coordinates
	r~{n}[] = Hypot[ X[]-cx~{i}~{j} , Y[]-cy~{i}~{j} ];
	n += 1;	
    EndFor 
EndFor

// Calculate MVP correction

For i In {1:NumWires}
     
    //Correction Terms for post-processing correction
    nu~{i}[] = (r~{i}[]>R_c) ? 1/(mur*mu0) :2* ((1/mu0)*(1 + mu0*((2*i[]*omega*sig*R_c^3*J_2[k[]*R_c])/(R_c*k[]^2*(J_0[k[]*R_c] - J_2[k[]*R_c])+2*k[]*J_1[k[]*R_c]))*A_s))/A_c;
   
EndFor
   
  dir[] =  sourceType==0 ? 1:-1;
  
}


//--------------------------------------------------------------------------------------------------------------
// Jacobian: Note that the Jacobian is set up for Full Model Integrals
//--------------------------------------------------------------------------------------------------------------

Jacobian {
  { Name Vol ;
    Case {
       { Region LR_DomainC; Jacobian Lin; Coefficient CoefJac[];}
       
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
	      { Region LR_Anode~{i};    Value Iin~{i}[]; }
          EndFor
      EndIf
    }
  }

  { Name Voltage_2D ;
    Case {
      If(!isSeries && sourceType == 1)
          For i In {1:NumWires}
	      { Region LR_Anode~{i};    Value Vin~{i}[]; }
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
	{ Region VoltageSource1; Value -V0; }
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
	     { Region LR_Anode~{i};     Branch{node_1, 1}; }
	  Else
	    { Region LR_Anode~{i};     Branch{node_1, node_2}; } 	
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
        In DomainCC; Jacobian Vol; Integration I1; }
      Galerkin { [nu_effective[] * Dof{d a} , {d a} ];
        In voidRegion; Jacobian Vol; Integration I1; }
          

      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }
      Galerkin { [ sigma[] * Dof{v} , {a} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }
      Galerkin { DtDof[sigma[] * Dof{a} , {v} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }


	 Galerkin { [ sigma[] * Dof{v} , {v} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }

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
      SetGlobalSolverOptions["-petsc_prealloc 3000"];
            
      Generate[S]; 
      Solve[S]; 
      
      // Postprocessing to take proximity losses into account
      PostOperation[SleeveArea];
      PostOperation[CoilCurrent];
      
      PostOperation[TruncatedFlux];
      PostOperation[FluxCorrection];
      PostOperation[Flux];
      
      PostOperation[Bav];
      PostOperation[ProximityLoss];
      
      SaveSolution[S]; 
      
      For i In {1:NumWires}
	Print[ {i, $sleeve_area~{i}, Sqrt[$sleeve_area~{i}/Pi], $fluxcorr~{i}, $flux~{i}, $bav~{i}^2*p_prox[], $P_prox~{i}, $coilcurrent~{i}},  
	       Format " Turn %g: sleeveArea = %7.5e, sleeveRadius = %7.5e, fluxcorr = %7.5e, flux = %7.5e, bav^2 = %7.5e , prox_loss = %7.5e, coilcurrent = %7.5e "] ;
      EndFor
      
       If(!isSeries)
	PostOperation[ProximityLossMatrix];
          For j In {1:NumWires}	  
              If(ActiveSource~{j} !=0 )		
                  For i In {1:NumWires}
	              If(i != j)
		          Print[ {i,j,$P_prox~{i}~{j} }, Format " P_prox_%2g_%2g = %7.3e"  ] ;
   	              EndIf
	          EndFor
              EndIf
          EndFor
      EndIf
      
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

{ Name az;  Value { Term { [ dir[]*CompZ[{a}]];     In Domain ; Jacobian Vol; } } }



//--------------------------------------------------------------------------------------------------------------
// Global Quantities:
//              reU: Real component of Voltage (V)
//              reI: Real component of Current (A)
//              imU: Imaginary component of Voltage (V)
//              imI: Imaginary component of Current (A)
//--------------------------------------------------------------------------------------------------------------
     
      
     // Retreives the resistance, self and mututual inductance in the case of parallel wires (inductance matrix)
     If(!isSeries)
         For j In {1:NumWires}	  
	     If(ActiveSource~{j} !=0 )		
	         For i In {1:NumWires}  
		     If(i==j)
		       
	                 { Name R~{i}~{j}; Value { Term { [ Re[{U}/{I}]];        In LR_Anode~{j}; } } }
		         { Name L~{i}~{j}; Value { Term { [ Im[{U}/{I}]/omega];  In LR_Anode~{j}; } } }
		    
	             Else
		      
		         { Name R~{i}~{j}; Value { Term { [ Re[i[]*omega*$flux~{i}/{I}]];       In LR_Anode~{j}; } } }
                         { Name L~{i}~{j}; Value { Term { [ Im[i[]*omega*$flux~{i}/{I}]/omega]; In LR_Anode~{j}; } } }                  
		      
		     EndIf
                 EndFor	    
	     EndIf
	 EndFor	  
     EndIf  
  
    
     // retreive total current/voltage in the circuit
     If(isSeries)
       
       { Name Rz;  Value { 
	   Term { [ Re[{Uz}/{Iz}] ];  In Sources_Cir; }
	   
	    
	   
	    
	   For i In {1:NumWires}
	     Term { [ Re[($P_prox~{i})/({Iz}*Conj[{Iz}])] ];  In Sources_Cir; }
	   EndFor
	   
	   } 
       }
       
       { Name Lz;  Value { 
	   
	   //Term { [ Im[{Uz}/{Iz}]/omega ];  In Sources_Cir; }      
	    
	   
	   For i In {1:NumWires}
	     Term { [ (Fabs[Re[$flux~{i}/{Iz}]]) ];   In Sources_Cir; }
	   EndFor
	   
	   
	   
	   } 
       }
       
       //----
       
       
         { Name reUz;  Value { 
	     Term { [ Re[{Uz}] ];  In Sources_Cir; } 
	   } 
	 }
	 //----
	 
         { Name reIz;  Value { 
	     Term { [ Re[{Iz}] ];  In Sources_Cir; } 
	   }
	 }
	 
         { Name imUz;  Value { 
	       Term { [ Im[{Uz}] ]; In Sources_Cir; }
	   } 
	 }
	 
         { Name imIz;  Value { 
	     Term { [ Im[{Iz}] ];  In Sources_Cir; } 
	   } 
	 }
	 
     EndIf

     

//-------------------------
    
    // Take proximity losses into account 
    //Area, flux correction, Average Flux density, etc
    For i In {1:NumWires} 
           { Name area~{i}; Value{ Integral{ [ 1 ] ;
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getTruncatedFlux~{i}; Value { Integral { [ CompZ[ {a} ]/$sleeve_area~{i}  ]; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getFluxCorrection~{i};  Value { 
	       Integral { [ (mu0 *  $coilcurrent~{i} / (2*Pi) * ( i[] * (skin_depth/R_c)^2 + Log[(Sqrt[$sleeve_area~{i}/Pi])/R_c] - J_0[k[]*R_c] / J_1[k[]*R_c] / (k[]*R_c) )) /$sleeve_area~{i}  ];  
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getFlux~{i}; Value { Integral { [ ($truncflux~{i}+$fluxcorr~{i})/$sleeve_area~{i} ]; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	 { Name getBav~{i}; Value {  Integral { [  (Norm[{d a}])/$sleeve_area~{i} ]; 
	       In  Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}} 
	 
	   { Name getProxLoss~{i}; Value { Integral {[ ($bav~{i}^2*p_prox[])/$sleeve_area~{i}   ] ; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; } } }
	 
	   { Name getCoilCurrent~{i};  Value { Term { [{I}];  
	       In LR_DomainC~{i}; Integration I1 ; Jacobian Vol; } } }
	 
    EndFor

    If(isSeries)
        { Name getCircuitCurrent;  Value { Term { [ {Iz} ];  
	       In Sources_Cir; } } }
    EndIf

     
     
     //--------------------------------------------------------------------------------------
     
         If(!isSeries)
    For j In {1:NumWires}	  
	    If(ActiveSource~{j} !=0 )		
	        For i In {1:NumWires}
		  If(i != j)
		    { Name getProxLossMatrix~{i}~{j};
                        Value { Integral {[ ($bav~{i}/$sleeve_area~{i})^2*p_prox[]/$sleeve_area~{i}   ] ; 
	                In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }	  	
                      }
                    } 
		  EndIf
		EndFor
            EndIf
	 EndFor
    EndIf
    
    
           If(!isSeries) 
       	// PROXIMITY EFFECT
       For j In {1:NumWires}	  
      If(ActiveSource~{j} !=0 )		
          For i In {1:NumWires}
	      If(i != j)
		{ Name P_proximity~{i}~{j}; Value {
		    Term { [ $P_prox~{i}~{j} + 0*{I} ];               In LR_DomainC~{j}; } } }
		
		{ Name R_proximity~{i}~{j}; Value {
		    Term { [ $P_prox~{i}~{j}/({I}*Conj[{I}]) ];       In LR_DomainC~{j}; } } }
	      	
   	      EndIf
	   EndFor
       EndIf
    EndFor
    EndIf


     
     }
   } 
 }



//--------------------------------------------------------------------------------------------------------------
// PostOperation
//--------------------------------------------------------------------------------------------------------------

// To add name to the correction (structured or unstructured sleeve)
If(withRing)
    ext_sleeve = Sprintf["-structured"];
Else
    ext_sleeve = Sprintf["-unstructured"]; 
EndIf

PostOperation PostOperations UsingPost PostProcessings {
  If(showMaps)

	     
    Print [az   , OnElementsOf Domain,   File  StrCat[ResultsDir,"az_MagDyn.pos"] ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;",
		 "View[l].IntervalsType = 3;",
		 "View[l].NbIso = 30;",
		 "View[l].RaiseZ = 15000;",
	         "View[l].RangeType = 3;"],
	       File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly] ;
       
       Print [ az, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
	  Format Gmsh, File StrCat[ResultsDir,"Cut_az.pos"] ];
          Echo[ StrCat["l=PostProcessing.NbViews-1;",
	       "View[l].Axes = 3;",
	       "View[l].LineWidth = 3;",
	       "View[l].Type = 2;"],
	File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly];
     
  EndIf

        Print [ az, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
		Format SimpleTable, File StrCat[ResultsDir,"Cut_az_thin",ext_sleeve,".dat"] ];
	      

  //Print[ flux_LMcorr[ LR_DomainC ], OnGlobal,
  //	 Format Table, File StrCat[ResultsDir,"flux_LMcorr.txt"], SendToServer StrCat["Output/0Impedance/1Inductance at Source (2D)/3Mat-FLUX"]];
       
       
	// Total voltage and current measured from source
	If(isSeries)
	  
	    Print [Rz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Rz.dat"],
	      SendToServer StrCat["Output/0Impedance/0Resistance at Source (2D)/3Mat",ext_sleeve], Color "AliceBlue"];
	    Print [Lz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Lz.dat"] ,
	      SendToServer StrCat["Output/0Impedance/1Inductance at Source (2D)/3Mat",ext_sleeve], Color "AliceBlue"];
	    
	    
	  
     Else
	 
		      // Retreive resistance, self and mutual inductance
	      For j In {1:NumWires}
                  If(ActiveSource~{j} !=0 )
	              For i In {1:NumWires}
		       
	                  ext_ij = StrCat[Sprintf["%g",i],Sprintf["%g",j]]; // Zji 
	                  If(i==j)
	                      Print[ R~{i}~{j}, OnRegion  LR_Anode~{j}, File > StrCat[ResultsDir,"Rper_full.dat"], Format Table,
		                  SendToServer StrCat[outputname_resistance_wire_thin,ext_sleeve, ext_ij] ];		   

	                  EndIf
	       
	                      Print[ L~{i}~{j}, OnRegion  LR_Anode~{j}, File > StrCat[ResultsDir,"Lper_full.dat"], Format Table,
		                  SendToServer StrCat[outputname_inductance_wire_thin,ext_sleeve, ext_ij] ];
		       
	              EndFor
	          EndIf
	      EndFor
	 
       EndIf
	
	

}



// Post-processing to capture proximity losses

PostOperation SleeveArea UsingPost PostProcessings {
    For i In {1:NumWires}
        Print[ area~{i}[Dom_Sleeve~{i}], OnGlobal, Format Table, File "", StoreInVariable $sleeve_area~{i}];
    EndFor 
}

PostOperation CoilCurrent UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getCoilCurrent~{i}, OnRegion LR_DomainC~{i}, 
       Format Table, File "", StoreInVariable $coilcurrent~{i}];
   EndFor
}


PostOperation TruncatedFlux UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getTruncatedFlux~{i}[Dom_Sleeve~{i}], OnGlobal, 
       Format Table, File "", StoreInVariable $truncflux~{i}];
   EndFor
}


PostOperation FluxCorrection UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getFluxCorrection~{i}[Dom_Sleeve~{i}], OnGlobal,
       Format Table, File "", StoreInVariable $fluxcorr~{i}];
   EndFor
}


PostOperation Flux UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getFlux~{i}[Dom_Sleeve~{i}], OnGlobal,
       Format Table, File "", StoreInVariable $flux~{i}];
   EndFor
}

PostOperation Bav UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getBav~{i}[Dom_Sleeve~{i}], OnGlobal,
       Format Table, File "", StoreInVariable $bav~{i}];
   EndFor
}

PostOperation ProximityLoss UsingPost PostProcessings {
   For i In {1:NumWires}
     Print[ getProxLoss~{i}[Dom_Sleeve~{i}], OnGlobal, 
       Format Table, File "", StoreInVariable $P_prox~{i}];
   EndFor
}


 If(!isSeries)
PostOperation ProximityLossMatrix UsingPost PostProcessings {
  For j In {1:NumWires}	  
      If(ActiveSource~{j} !=0 )		
          For i In {1:NumWires}
	      If(i != j)
                  Print[ getProxLossMatrix~{i}~{j}[Dom_Sleeve~{i}], OnGlobal, 
                      Format Table, File "", StoreInVariable $P_prox~{i}~{j}];
   	      EndIf
	   EndFor
       EndIf
  EndFor
}
EndIf





DefineConstant[
  R_ = {"Resolutions"    , Name "GetDP/1ResolutionChoices"      , Visible 1},
  C_ = {"-solve -v2"     , Name "GetDP/9ComputeCommand"         , Visible 1},
  P_ = {"" , Name "GetDP/2PostOperationChoices"   , Visible 1}
];
