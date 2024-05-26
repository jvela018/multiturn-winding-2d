Include "Common/parameters.pro";


// Initial parameters and point mesh control 
If(isThin)
    rad = rs;
Else
    rad = R_c;
EndIf
If(isThin)
    msh_lconductor = rad;
    msh_voutline = rad/3;
Else
    If(freqRefinement)
        msh_lconductor = mesh_parameter;
	msh_voutline = mesh_parameter;
    Else
        msh_lconductor = rad/10;
	msh_voutline = rad/3;
    EndIf
EndIf

msh_air = 3e-3;


// Create Points as Thin Wires
// Define Points For Round/Squared Outline (withRing) 
For j In {1:columns} 
    For i In {1:rows}
          
        // Center points are the thin wires on physical entities 
        CenterPoints[] += newp; Point(newp) = {cx~{i}~{j}, cy~{i}~{j}, 0, msh_lconductor};	  
    
        If(withRing)
    
            If(isRound)
                //Points that create the round conductor outline---------------------------------------------------------
		OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}+rad, cy~{i}~{j}, 0, msh_voutline};	       
	        OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}, cy~{i}~{j}+rad, 0, msh_voutline};
		
		// If conductors are touching in -x, create common point else create separate points for each conductor
		If(j == 1 || distx > 2*rad)
		    OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}-rad, cy~{i}~{j}, 0, msh_voutline};
		EndIf
	       
	      
		// If conductors are touching in -y, create common point else create separate points for each conductor
		If(i == 1 || disty > 2*rad)
		    OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}, cy~{i}~{j}-rad, 0, msh_voutline};
		EndIf
		
		
	        //-------------------------------------------------------------------------------------------------------
            Else
	        //Points that create the squared conductor outline ------------------------------------------------------
	        OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}+rad, cy~{i}~{j}+rad, 0, msh_voutline};
	        OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}-rad, cy~{i}~{j}+rad, 0, msh_voutline};
		
		// If conductors are touching in -x, create common point else create separate points for each conductor
		If(j == 1 || distx > 2*rad)
	            OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}-rad, cy~{i}~{j}-rad, 0, msh_voutline};
		EndIf
		
		// If conductors are touching in -y, create common point else create separate points for each conductor
		If(i == 1 || disty > 2*rad)
	            OutlinePoints[] += newp; Point(newp) = {cx~{i}~{j}+rad, cy~{i}~{j}-rad, 0, msh_voutline};
		EndIf
	        //-------------------------------------------------------------------------------------------------------
            EndIf
    
        EndIf
		 
    EndFor
EndFor    


//-----------------------------------------------------------------------------------------------------------------------
// Create Lines For Conductor Outline
//-----------------------------------------------------------------------------------------------------------------------

If(withRing)
    line_index = 0;
    linet_index = 0;
   
    For i In {0:#CenterPoints[]-1} 
       For j In {0:3} //   For every center point 4 lines are defined
      
	   If(isRound) //Round Conductor Outline
	       	      
               If(j < 3)
                   OutlinesLines[] += newl; Circle(newl) = {OutlinePoints[j+line_index], CenterPoints[i], OutlinePoints[j+1+line_index]};
               Else
	           OutlinesLines[] += newl; Circle(newl) = {OutlinePoints[j+line_index], CenterPoints[i], OutlinePoints[0+line_index]};
               EndIf
	
	   Else // Squared Conductor Outline
	       If (j <3 )
                    OutlinesLines[] += newl; Line(newl) = {OutlinePoints[j+line_index], OutlinePoints[j+1+line_index]};	
               Else
                    OutlinesLines[] += newl; Line(newl) = {OutlinePoints[j+line_index], OutlinePoints[0+line_index]};
               EndIf  	 
	   EndIf	 
	   
       EndFor
     
    
    // Create conductor surfaces
    ConductorLineLoop[] += newll; Curve Loop(newll) = {OutlinesLines[{0+line_index,1+line_index,2+line_index,3+line_index}]};
    ConductorSurface[] += news; Plane Surface(news) = {ConductorLineLoop[i]};
    line_index += 4;
    
    //Embedd Center-line Points in Mesh
    Point{CenterPoints[i]} In Surface{ConductorSurface[i]};
    
    EndFor
EndIf


If(AirAxi)
    AirPoints[] += newp; Point(newp) = {0, -air_distance, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {0, 0, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {0, air_distance, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {air_distance, 0, 0, msh_air};

    AirLines[] += newl; Line(newl)   = {AirPoints[0], AirPoints[1]};
    AirLines[] += newl; Line(newl)   = {AirPoints[1], AirPoints[2]};
    AirLines[] += newl; Circle(newl) = {AirPoints[2], AirPoints[1], AirPoints[3]};
    AirLines[] += newl; Circle(newl) = {AirPoints[3], AirPoints[1], AirPoints[0]};
Else
      
    AirPoints[] += newp; Point(newp) = {-0.1, -0.1, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {0.1, -0.1, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {0.1, 0.1, 0, msh_air};
    AirPoints[] += newp; Point(newp) = {-0.1, 0.1, 0, msh_air};

    AirLines[] += newl; Line(newl)   = {AirPoints[0], AirPoints[1]};
    AirLines[] += newl; Line(newl)   = {AirPoints[1], AirPoints[2]};
    AirLines[] += newl; Line(newl)   = {AirPoints[2], AirPoints[3]};
    AirLines[] += newl; Line(newl)   = {AirPoints[3], AirPoints[0]};
   
      
EndIf

// Create conductor surfaces
AirLineLoop[] += newll; Curve Loop(newll) = {AirLines[{0,1,2,3}]};

If(withRing)
    AirSurface[] += news; Plane Surface(news) = {AirLineLoop[0],ConductorLineLoop[]};
Else
    AirSurface[] += news; Plane Surface(news) = {AirLineLoop[0]};
    
    //Embedd Center-line Points in Mesh
    For i In {0:#CenterPoints[]-1} 
        Point{CenterPoints[i]} In Surface{AirSurface[0]};
    EndFor
EndIf

// ------------------------------------------------------------------------
// Mesh controls: these re-write initial point mesh control
// ------------------------------------------------------------------------

// If is thin then outline becomes sleeve domain
If(isThin && withRing)
    sleeve_mesh_index = 0 ;
    For i In {1:#CenterPoints[]}  
      Transfinite Curve {OutlinesLines[{0+sleeve_mesh_index,1+sleeve_mesh_index,2+sleeve_mesh_index,3+sleeve_mesh_index}]} = 3 Using Progression 1;//135
        sleeve_mesh_index += 4;
    EndFor
EndIf

If(!isThin  && !freqRefinement)
    sleeve_mesh_index = 0 ;
    For i In {1:#CenterPoints[]}  
      Transfinite Curve {OutlinesLines[{0+sleeve_mesh_index,1+sleeve_mesh_index,2+sleeve_mesh_index,3+sleeve_mesh_index}]} = 80 Using Progression 1; //40 135 and 100 for far distance?
        sleeve_mesh_index += 4;
    EndFor
EndIf

// Define boundary layers / mesh control for volume conductors
If(!isThin && isBoundaryLayer && !freqRefinement)   
    Field[1] = BoundaryLayer;
    Field[1].EdgesList =  {Boundary{ Surface{ConductorSurface[]};}};
    Field[1].hfar = rad/8;                      // Element size far from the wall
    Field[1].hwall_n = rad/5000;                 // Mesh Size Normal to the The Wall
    Field[1].thickness = rad/7;                  // Maximum thickness of the boundary layer
    Field[1].ratio = 1.07;                       // Size Ratio Between Two Successive Layers
    Field[1].ExcludedFaceList = {AirSurface[]};  // Don't create BL in this surface
    Field[1].Quads = 1;                          // Make quads
    BoundaryLayer Field = 1;
EndIf

// Air mesh control
If(AirAxi)
  Transfinite Curve {AirLines[{0,1}]} = 15 Using Progression 1;
  Transfinite Curve {AirLines[{2,3}]} = 15 Using Progression 1;
Else
    Transfinite Curve {AirLines[{0,1,2,3}]} = 20 Using Progression 1; //Originally 50
EndIf


// ------------------------------------------------------------------------
// Define Physical Entities
// ------------------------------------------------------------------------

Physical Curve("Infinity", Inf) = {AirLines[{2,3}]};
Physical Curve("No Flux Symmetry Boundary", noFlux) = {AirLines[{0,1}]};
Physical Surface("Full Model - DomainCC ", Full_DomainCC) = {AirSurface[]};

// If is the full model create conductor surfaces
If(!isThin)
    For i In {0:#CenterPoints[]-1}
        Physical Surface(Sprintf("Full Model - DomainC_%g",i), Full_DomainC+i) = {ConductorSurface[i]};
    EndFor
EndIf

// If it's semi-analytical model create physical entities depending on structure/unstructured sleeve

// structured sleeve belongs to DomainCC
If(isThin && withRing) 
    If(isThin == 3)
        For i In {0:#CenterPoints[]-1}
	  Physical Surface(Sprintf("Full Model - DomainC_%g",i), Full_DomainC+i) = {ConductorSurface[i]};
	EndFor
    Else
        Physical Surface("Full Model - DomainCC ", Full_DomainCC) += {ConductorSurface[]};
    EndIf
EndIf

// Create thin wires
If(isThin)
    For i In {0:#CenterPoints[]-1}
        Physical Point(Sprintf("Thin Model - DomainC_%g",i), LR_DomainC+i) = {CenterPoints[i]};
    EndFor
EndIf

