//--------------------------------------------------------------------------------------------------------------
// Full/Thin 2D Model/Mesh
// EM Formulations: Magnetodynamics
// Multi-turn coil Model
//--------------------------------------------------------------------------------------------------------------

Include "Common/parameters.pro";

// 2D Models

If(isThin == 0)
    Include "Models/full_2D.pro";
ElseIf(isThin==2)
    Include "Models/thin_semi_analytical_2D.pro";
ElseIf(isThin==1)
    Include "Models/thin_no_correction_2D.pro";
    //ElseIf(isThin==3)
    //Include "Models/thin_material_2D.pro";
EndIf



