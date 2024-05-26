//--------------------------------------------------------------------------------------------------------------
// Full 2D Single Conductor Model/Mesh
// EM Formulations: Magnetodynamics
// Simple example using Global Sources (Voltage/Current)
// with and without circuit connections
//--------------------------------------------------------------------------------------------------------------


// Create Folders for Results
 DefineConstant[FileId = ""];
DefineConstant[ResultsDir  = StrCat["getdp_results_thin_2D/",FileId]];


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
	Full_DomainC  += Region[{Full_DomainC~{i+1}}];
    EndFor

  
  If(withIED)
      IED = Region[{IED}];
      DomainCC += Region[{IED,Full_DomainCC}];  
  Else
      DomainCC += Region[{Full_DomainCC, Full_DomainC}];  
  EndIf
  
  
  Domain   = Region[{DomainC, DomainCC}];
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
      Dom_Sleeve~{i} = ElementsOf[ DomainCC, OnOneSideOf LR_DomainC~{i} ];
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

  sigma[LR_DomainC]   = sig;      //Electrical Conductivity
  sigma[DomainCC]  = 0.;       //Electrical Conductivity in Non-Conducting Domain

  mu[LR_DomainC]      = mur*mu0;  //Permeability of Free Space
  mu[DomainCC]     = mur*mu0;  //Permeability of Conducting Domain
  nu[]             = 1/mu[];   //Reluctivity

//--------------------------------------------------------------------------------------------------------------
// Excitation Parameters
//--------------------------------------------------------------------------------------------------------------
  i[] = Complex[0,1];
  Iin[] = I0;
  Vin[] = V0;

  //Actual Area of the conductors/Sleeve
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
  tau[] = (i[]-1)*R_c/skin_depth;


  L_int[] = Im[(R_dc/2)*k[]*R_c*J_0[k[]*R_c]/J_1[k[]*R_c]]/omega; // 
  R_skin[] = Re[(R_dc/2)*(k[]*R_c)*J_0[k[]*R_c]/J_1[k[]*R_c]];    // 
  Correction_flux[]   = mu0 * Iin[] / (2*Pi) * ( i[] * (skin_depth/R_c)^2 + Log[rs/R_c] - J_0[k[]*R_c] / J_1[R_c*k[]] /(k[]*R_c) ) ; 

 // Symmetry
 //SymmetryFactor = 2; //Note we are modeling only half of the circle
 

  //Prox stuff from Joule Loss
  p_denom[] = (R_c*k_conj[]*(J_0[k_conj[]*R_c] - J_2[k_conj[]*R_c]) + 2*J_1[k_conj[]*R_c])*(R_c*k[]*(J_0[k[]*R_c] - J_2[k[]*R_c]) + 2*J_1[k[]*R_c]) * (k[]^2 - k_conj[]^2);
  p_nume[] = (16*omega^2*sig*R_c^2)*(Pi*R_c*(k_conj[]*J_1[k[]*R_c]*J_0[k_conj[]*R_c] - k[]*J_0[k[]*R_c]*J_1[k_conj[]*R_c]));
  p_prox[] = p_nume[]/p_denom[];
 
 
  J[] = (Iin[]/A_c)*Vector[0,0,1];
  dir[] =  sourceType==0 ? 1:-1;
  // Derivative of J1 evaluated at Rc
  dJ_1_R[] = 0.5*k[]*(J_0[k[]*R_c] - J_2[k[]*R_c]);
  
  CoefJac[]=A_c;
  
  DefineConstant[ skind   = {skin_depth, Name "Input/Skin Effect/0Skin Depth (m) " , Highlight "LightGreen", Visible 1, ReadOnly}
                  condrad = {R_c  , Name "Input/Skin Effect/1Conductor Radius (m) ", Highlight "LightGreen", Visible 1, ReadOnly}
		  skin_rad_ratio = {R_c/skin_depth  , Name "Input/Skin Effect/2Ratio", Highlight "LightGreen", Visible 1, ReadOnly}
		  mesh_para = {mesh_parameter  , Name "Input/Skin Effect/3 Mesh Characteristic Length", Highlight "LightGreen", Visible 0, ReadOnly}
  ];
  
  // set up voltage/current input sources depending on parallel or series connection
    For i In {1:NumWires}
      If(isSeries)
          Vin~{i}[] = I0;
          Iin~{i}[] = V0;
      Else
	  Vin~{i}[] = Vp~{i};
          Iin~{i}[] = Ip~{i};  
      EndIf
  EndFor
  
// cylindrical r coordinate for correction of magnetic vector potential in post-processing 
n = 1;
For i In {1:rows}
    For j In {1:columns}
        // converting x and y coordinates to r coordinates
	r~{n}[] = Hypot[ X[]-cx~{i}~{j} , Y[]-cy~{i}~{j} ];
	theta~{n}[] = Atan2[ Y[]-cy~{i}~{j} , X[]-cx~{i}~{j} ];
	n += 1;	
    EndFor 
EndFor

// Calculate MVP correction
B0 = 0.1;
For i In {1:NumWires}
  
  //Analytical_az_Prox~{i}[] =  -1e12 * Cos[theta~{i}[]] * ( (r~{i}[]>R_c) ? 
  //    R_c^2/r~{i}[] * ( 2*J_1[ tau[] ] / ( tau[] * J_0[ tau[] ] ) - 1 ) :
  //    ( 2*R_c * J_1[ tau[]*r~{i}[]/R_c ] / ( tau[] * J_0[ tau[] ] ) - r~{i}[] ) );
    
  //Analytical_az_Prox~{i}[]  = (r~{i}[] <= R_c) ? ((-2*Norm[$by~{i}^2 +$bx~{i}^2] *R_c)/(R_c*dJ_1_R[] + J_1[k[]*R_c]))*J_1[k[]*r~{i}[]]*(Sin[theta~{i}[]]) : 0; 
  
   Analytical_az_Prox~{i}[] =  -$by~{i} * Cos[theta~{i}[]] * ( (r~{i}[]>R_c) ? 
      R_c^2/r~{i}[] * ( 2*J_1[ tau[] ] / ( tau[] * J_0[ tau[] ] ) - 1 ) :
      ( 2*R_c * J_1[ tau[]*r~{i}[]/R_c ] / ( tau[] * J_0[ tau[] ] ) - r~{i}[] ) );
    
    // analytical corrections (wrt the truncated field) 
    Analytical_az_Skin~{i}[]  = mu0 *I0 / (2*Pi) * ( (r~{i}[]>R_c) ? Log[rs/r~{i}[]] :Log[rs/R_c] + (J_0[k[]*r~{i}[]] - J_0[k[]*R_c]) / J_1[k[]*R_c] / (k[]*R_c) );
      
    Analytical_a~{i}[] = Analytical_az_Skin~{i}[] + Analytical_az_Prox~{i}[];
    //Correction Terms for post-processing correction
    Correction_a~{i}[]    = ((r~{i}[]< rs) ? Analytical_a~{i}[] : 0)* Vector[0,0,1];
   
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
	{ Region CurrentSource1; Value -Iin~{i}[]; }
      EndIf
    }
  }
  { Name Voltage_Cir ;
    Case {
      If(isSeries && sourceType == 1)
	{ Region VoltageSource1; Value Vin~{i}[]; }
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
  // 2D

  { Name Hcurl_a_2D; Type Form1P;
    BasisFunction {
      { Name sc; NameOfCoef ac; Function BF_PerpendicularEdge;
	Support Dom_Hcurl_a; Entity NodesOf[All]; }
      { Name sw; NameOfCoef aw; Function BF_PerpendicularEdge;
	Support Dom_Hthin_a; Entity NodesOf[LR_DomainC]; }
    }
    SubSpace {
      { Name ac ; NameOfBasisFunction { sc } ; }
      { Name aw ; NameOfBasisFunction { sw } ; }
    }
    Constraint {
      { NameOfCoef ac;
	EntityType Auto; NameOfConstraint MVP_2D; }
    }
  }


  
  // Current in stranded coil (2D)
  { Name Hregion_i_2D; Type Vector;
    BasisFunction {
      { Name sr; NameOfCoef ir; Function BF_RegionZ;
        Support LR_DomainC; Entity LR_DomainC; }
    }
    GlobalQuantity {
      { Name I; Type AliasOf; NameOfCoef ir; }
      { Name U; Type AssociatedWith; NameOfCoef ir; }
    }
    Constraint {
      { NameOfCoef U;
        EntityType Region; NameOfConstraint Voltage_2D; }
	
      { NameOfCoef I;
        EntityType Region; NameOfConstraint Current_2D; }
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
      { Name a;  Type Local;  NameOfSpace Hcurl_a_2D; }
      { Name ac; Type Local;  NameOfSpace Hcurl_a_2D[ac]; }
      { Name aw; Type Local;  NameOfSpace Hcurl_a_2D[aw]; }
      
      { Name i; Type Local;  NameOfSpace Hregion_i_2D; }
      { Name U;  Type Global; NameOfSpace Hregion_i_2D [U]; }
      { Name I;  Type Global; NameOfSpace Hregion_i_2D [I]; }
      
      If(isSeries)
          { Name Uz; Type Global;   NameOfSpace Hregion_Cir [Uz]; }
          { Name Iz; Type Global;   NameOfSpace Hregion_Cir [Iz]; }
      EndIf
      
    }

    Equation {

      
      Integral { [ nu[]* Dof{d ac} , {d ac} ];
        In Domain; Jacobian Vol; Integration I1; }    
      //Integral { [  -J[] , {ac} ];
      //In LR_DomainC; Jacobian Vol; Integration I1; }
      
      Integral { [ nu[] * Dof{d aw} , {d aw} ];
        In Domain; Jacobian Vol; Integration I1; }
      //Integral { [ -J[] , {aw} ];
      //In LR_DomainC; Jacobian Vol; Integration I1; }
      
    
      GlobalTerm { [-Dof{U} , {I} ]; In LR_DomainC; }    
    
    
      For i In {1:NumWires}
	GlobalTerm { [ R_skin[] * Dof{I} , {I} ]; 
	  In LR_DomainC~{i};} 
	GlobalTerm { DtDof [ L_int[] * Dof{I} , {I} ]; 
	In LR_DomainC~{i};} 
       EndFor
      
       //Coupling
      Integral {  [ -Dof{i}/A_c , {ac} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }   
      Integral {  [ -Dof{i}/A_c , {aw} ];
        In LR_DomainC; Jacobian Vol; Integration I1; }
	 	        
      Integral { DtDof [ Dof{ac}/A_c , {i} ];
       In LR_DomainC; Jacobian Vol; Integration I1; }   
      Integral { DtDof [ -Dof{aw}/A_c , {i} ];
       In LR_DomainC; Jacobian Vol; Integration I1; }
	

     If(isSeries)
         GlobalTerm{ [ 0 * Dof{Iz},       {Iz} ];  In Domain_Cir; }

         // Inserting the network
         GlobalEquation { Type Network; NameOfConstraint ElectricalCircuit;
             { Node {I};     Loop {U};     Equation {U};   In Electrodes_with_GQ; }
             { Node {Iz};    Loop {Uz};    Equation {Uz};  In Domain_Cir; }
          }
     EndIf
      
      
      
      //---------------------------------------------------------------------
      
      
      //---------------------------------------------------------------------

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
      
      // Print number of Degrees-of-Freedom onto File
      Print[{Freq, $KSPSystemSize, nElements[]}, Format "%g %g %g", File StrCat[ResultsDir, "2D_Dof_nElem_Correction.dat"] ];
      
      // Postprocessing to take proximity losses into account
      PostOperation[SleeveArea];
      PostOperation[CoilCurrent];
      
      PostOperation[TruncatedFlux];
      PostOperation[FluxCorrection];
      PostOperation[Flux];
      
      PostOperation[Bav];
      PostOperation[ProximityLoss];
      
      PostOperation[b_proximity];
      For i In {1:NumWires}
      Evaluate[ $bx~{i} = $bxInteg~{i} / $sArea~{i} ];
      Evaluate[ $by~{i} = $byInteg~{i} / $sArea~{i} ];
      Print[ { i, $sArea~{i}, $bx~{i}, $by~{i} },
          Format "Conductor %g: area=%7.5e bx=%7.5e by=%7.5e"];
      EndFor
      
      
      
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

{ Name acwz;
	Value { Local { [ CompZ[ {ac} - {aw} ] ];
	    In  Dom_Hcurl_a; Jacobian Vol; }}}

{ Name ac;  
  Value { 
    Term { [ (CompZ[{ac}]) ];In Domain ; Jacobian Vol; } 

  } 
}

{ Name az; Value {
    Local { [ dir[]*(CompZ[{ac}] - CompZ[{aw}]) ]; In  Dom_Hcurl_a; Jacobian Vol; } 
        
    For i In { 1:NumWires}
      Term {[CompZ[dir[]*Correction_a~{i}[] ]];In Dom_Hthin_a; Jacobian Vol; } 
    EndFor
  } 
}

{ Name jz;  Value {
    Term { [ -J[] ]; In DomainC ; Jacobian Vol; } 
  } 
}

{ Name flux_LMcorr;
  Value { 
             For i In {1:NumWires}
	       Integral { [ (CompZ[ {ac} - {aw} ] + Correction_flux[])/(Iin[]) ]; 
		     In LR_DomainC~{i}; Integration I1 ; Jacobian Vol; }
	     EndFor
   }
}

{ Name acorr_skin;  
  Value { 
    For i In { 1:NumWires}
        Term { [ ( Correction_a~{i}[] ) ];In Dom_Hcurl_a; Jacobian Vol; } 
    EndFor
  } 
}

{ Name acorr_prox;  
  Value { 
    For i In { 1:NumWires}
      Term { [ (Correction_a~{i}[] ) ];In Dom_Hcurl_a ; Jacobian Vol; } 
    EndFor
  } 
}


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
			 { Name L~{i}~{j}; Value { Term { [ Im[i[]*$flux~{i}/{I}]]; In LR_Anode~{j}; } } }                  
		      
		     EndIf
                 EndFor	    
	     EndIf
	 EndFor	  
     EndIf  
  
    
     // retreive total current/voltage in the circuit
     If(isSeries)
       
       { Name Rz;  Value { 
	   Term { [ Re[-{Uz}/{Iz}] ];  In Sources_Cir; }
	   
	   For i In {1:NumWires}
	     Term { [ Re[($P_prox~{i})/({Iz}*Conj[{Iz}])]  ];  In Sources_Cir; }
	   EndFor
	   
	   
	   } 
       }
       
       { Name Lz;  Value { 
	   
	   //Term { [ -Im[{Uz}/{Iz}]/omega ];  In Sources_Cir; }      
	    
	   
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
	 
	   { Name getTruncatedFlux~{i}; Value { Integral { [ CompZ[ {ac} - {aw} ]/$sleeve_area~{i}  ]; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getFluxCorrection~{i};  Value { 
	       Integral { [ (mu0 *  $coilcurrent~{i} / (2*Pi) * ( i[] * (skin_depth/R_c)^2 + Log[(Sqrt[$sleeve_area~{i}/Pi])/R_c] - J_0[k[]*R_c] / J_1[k[]*R_c] / (k[]*R_c) )) /$sleeve_area~{i}  ];  
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getFlux~{i}; Value { Integral { [ ($truncflux~{i}+$fluxcorr~{i})/$sleeve_area~{i} ]; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}}
	 
	   { Name getBav~{i}; Value {  Integral { [  Norm[CompY[{d ac}-{d aw}]+CompX[{d ac}-{d aw}]]/$sleeve_area~{i} ]; 
	       In  Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; }}} 
	 
	   { Name getProxLoss~{i}; Value { Integral {[ ($bav~{i}^2*p_prox[])/$sleeve_area~{i}   ] ; 
	       In Dom_Sleeve~{i}; Integration I1 ; Jacobian Vol; } } }
	 
	   { Name getCoilCurrent~{i};  Value { Term { [{I}];  
	       In LR_DomainC~{i}; Integration I1 ; Jacobian Vol; } } }
	 
    EndFor
    
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

    
    For i In {1:NumWires}
        { Name SleeveArea~{i}; Value {
            Integral { [ 1 ];
              In Dom_Sleeve~{i}; Jacobian Vol; Integration I1;} } }
        { Name bx_integrated~{i}; Value {
            Integral { [ CompX[ {d a} ] ];
              In Dom_Sleeve~{i}; Jacobian Vol; Integration I1;} } }
        { Name by_integrated~{i}; Value {
            Integral { [ CompY[ {d a} ] ];
              In Dom_Sleeve~{i}; Jacobian Vol; Integration I1;} } }
      EndFor
     
     }
   } 
 }



//--------------------------------------------------------------------------------------------------------------
// PostOperation
//--------------------------------------------------------------------------------------------------------------
ExtGmsh    = ".pos" ;
// To add name to the correction (structured or unstructured sleeve)
If(withRing)
    ext_sleeve = Sprintf["-structured"];
Else
    ext_sleeve = Sprintf["-unstructured"]; 
EndIf

PostOperation PostOperations UsingPost PostProcessings {
  If(showMaps)
   Print [jz   , OnElementsOf DomainC,   File  StrCat[ResultsDir,"jz_MagDyn.pos"] ];

   Print [acwz   , OnElementsOf Domain,   File  StrCat[ResultsDir,"acwz_MagDyn.pos"] ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;",
		 "View[l].IntervalsType = 3;",
		 "View[l].NbIso = 30;",
		 "View[l].RaiseZ = 15000;",
	         "View[l].RangeType = 3;"],
	       File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly] ;
	     
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
      
      // Test for proximity test changing some stuff that no one cares about
       Print [ az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	  Format Gmsh, File StrCat[ResultsDir,"Cut_az.pos"] ];
          Echo[ StrCat["l=PostProcessing.NbViews-1;",
	       "View[l].Axes = 3;",
	       "View[l].LineWidth = 3;",
	       "View[l].Type = 2;"],
	File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly];
      
     
  EndIf

       Print [ az, OnLine {{-0.009,0,0}{0.009,0,0}}{500},
     Format Gmsh, File StrCat[ResultsDir, "a_line",ExtGmsh]];
   Echo[ StrCat["l=PostProcessing.NbViews-1;",
      "View[l].Axes = 3;",
      "View[l].LineWidth = 3;",
      "View[l].Type = 4;",
      "View[l].Name = Sprintf('a corrected (%g nodes)',Mesh.NbNodes);"], 
    File StrCat[ResultsDir,"tmp.pos"] ];
  
   Print [ acorr_prox, OnLine {{-0.009,0,0}{0.009,0,0}}{500},
     Format Gmsh, File StrCat[ResultsDir, "acorr_prox",ExtGmsh]];
   Echo[ StrCat["l=PostProcessing.NbViews-1;",
      "View[l].Axes = 3;",
      "View[l].LineWidth = 3;",
      "View[l].Type = 4;",
      "View[l].Name = Sprintf('a proxonly (%g nodes)',Mesh.NbNodes);"], 
    File StrCat[ResultsDir,"tmp.pos"] ];
  
  Print [ acorr_skin, OnLine {{-0.009,0,0}{0.009,0,0}}{500},
     Format Gmsh, File StrCat[ResultsDir, "acorr_skin",ExtGmsh]];
   Echo[ StrCat["l=PostProcessing.NbViews-1;",
      "View[l].Axes = 3;",
      "View[l].LineWidth = 3;",
      "View[l].Type = 4;",
      "View[l].Name = Sprintf('a skinonly (%g nodes)',Mesh.NbNodes);"], 
    File StrCat[ResultsDir,"tmp.pos"] ];
       
  Print [ acorr_prox, OnLine {{-0.009,0,0}{0.009,0,0}}{1000},
     Format SimpleTable, File StrCat[ResultsDir, "acorr_prox.txt"]];
   
   Print [ acorr_skin, OnLine {{-0.009,0,0}{0.009,0,0}}{1000},
     Format SimpleTable, File StrCat[ResultsDir, "acorr_skin.txt"]];
   
 
	      
       Print [ acwz, OnLine { {-0.009,0,0}{0.009,0,0} } {1000},
	    Format SimpleTable, File StrCat[ResultsDir,"Cut_acw.txt"] ];
	      
   Print [ ac, OnLine { {-0.009,0,0}{0.009,0,0} } {1000},
		Format SimpleTable, File StrCat[ResultsDir,"Cut_ac_thin",ext_sleeve,".txt"] ];
	      
   Print [ az, OnLine { {-0.009,0,0}{0.009,0,0} } {1000},
		Format SimpleTable, File StrCat[ResultsDir,"Cut_az_thin",ext_sleeve,".txt"] ];	      
  
  
    If(isSeries)
        Print [ az, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
		Format SimpleTable, File StrCat[ResultsDir,"Cut_az_thin",ext_sleeve,".dat"] ];
	      
        Print [ acwz, OnLine { {left_x,bottom_y,0} {right_x,top_y,0} } {800},
		Format SimpleTable, File StrCat[ResultsDir,"Cut_acwz_thin",ext_sleeve,".dat"] ];
  
	      //Print[ flux_LMcorr[ LR_DomainC ], OnGlobal,
	      // Format Table, File StrCat[ResultsDir,"flux_LMcorr.txt"], SendToServer StrCat["Output/0Impedance/1Inductance at Source (2D)/1Corr-FLUX"]];
     Else
	  Print [ az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	 Format SimpleTable, File StrCat[ResultsDir,"Cut_az_thin",ext_sleeve,".dat"] ];
     EndIf
	
		Echo[ StrCat["l=PostProcessing.NbViews-1;",
	       "View[l].Axes = 3;",
	       "View[l].LineWidth = 3;",
	       "View[l].Type = 2;"],
	       File StrCat[ResultsDir,"tmp.geo"], LastTimeStepOnly];
	     
	 	 Print [  az, OnLine { {-0.009,0,0} {0.009,0,0} } {800},
	  Format SimpleTable, File StrCat[ResultsDir,"Cut_az.dat"] ];
	      
       
       	// Total voltage and current measured from source
	If(isSeries)
	  
	    Print [Rz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Rz.dat"],
	      SendToServer StrCat["Output/0Impedance/0Resistance at Source (2D)/1Corr",ext_sleeve], Color "AliceBlue"];
	    Print [Lz,   OnRegion Sources_Cir,   Format Table, File > StrCat[ResultsDir,"Lz.dat"] ,
	      SendToServer StrCat["Output/0Impedance/1Inductance at Source (2D)/1Corr",ext_sleeve], Color "AliceBlue"];
	    
	    
	  
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
       
       
       If(!isSeries)
       // PROXIMITY EFFECT
       For j In {1:NumWires}	  
      If(ActiveSource~{j} !=0 )		
          For i In {1:NumWires}
	      If(i != j)
	          ext_ij = StrCat[Sprintf["%g",i],Sprintf["%g",j]]; // Zji 
		
		  //Print[ P_proximity~{i}~{j}, OnRegion  LR_DomainC~{j}, File > "P_prox_corr.dat", Format Table,
		  //   SendToServer StrCat[outputname_rePower_prox,"/1Correction-Prox", ext_ij] ];
		   
	        Print[ R_proximity~{i}~{j}, OnRegion   LR_DomainC~{j}, File > "R_prox_corr.dat", Format Table,
		  SendToServer StrCat[outputname_resistance_prox,ext_sleeve,"-Prox", ext_ij] ];
	      	
   	      EndIf
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


PostOperation b_proximity UsingPost PostProcessings {
  For i In {1:NumWires}
    Print[SleeveArea~{i}[ Dom_Sleeve~{i} ], OnGlobal, 
        Format TimeTable, File > "proximity.txt",
        StoreInVariable $sArea~{i}];
    Print[bx_integrated~{i}[ Dom_Sleeve~{i}], OnGlobal, 
      Format TimeTable, File > "proximity.txt",
      StoreInVariable $bxInteg~{i}];
    Print[by_integrated~{i}[ Dom_Sleeve~{i}], OnGlobal, 
      Format TimeTable, File > "proximity.txt",
      StoreInVariable $byInteg~{i}];
  EndFor
}



DefineConstant[
  R_ = {"Resolutions"    , Name "GetDP/1ResolutionChoices"      , Visible 1},
  C_ = {"-solve -v2"     , Name "GetDP/9ComputeCommand"         , Visible 1},
  P_ = {"" , Name "GetDP/2PostOperationChoices"   , Visible 1}
];
