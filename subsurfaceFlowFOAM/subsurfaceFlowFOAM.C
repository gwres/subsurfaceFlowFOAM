/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
/*-----------------------------------------------------------------------*\
Class
    Foam::subsurfaceFlowFOAM

Group
    subsurfaceFlowFOAM

Description
    Modeling of Subsurface Flow [Richards Equation]
\*-----------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    std::clock_t startT= std::clock(); 									//Start Time
    
    argList::addNote
    (
        "Subsurface Flow Model"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    
    //------------------------------------------------------------------------------------------------------------------------------//
    
    #include "readRichardParameters.H" 									//Reading the input parameters for "Richards Equation"
    
    #include "readPicardParameters.H" 									//Reading the input parameters for Picard Iteration Method
    
    //------------------------------------------------------------------------------------------------------------------------------//
        
    label N_ele = mesh.cells().size(); 									//Number of Elements
    
    label nSs = 0;            	       									//Indicator for zero Ss-field
    
    label nMARK = 0;           	       									//Indicator for Convergence Failure
    
    scalar errorSum = 0.0;       	   									//Initialization of Picard iteration error
    scalar error = 0.0;
    
    label sc = 0;				       									//Initialization of the Stabilization counter
        
    label NP = NPl;           	       									//Initialization of Picard Iteration number
    
    scalar MB2 = 0.0;			       									//Initialization of total net flux variable
    
    scalar oldTime = 0.0;	  	       									//Initialization of previous time variable
	label oldTimeIndex = 0;	  	       									//Initialization of previous time-index variable
    
    volScalarField h0 = h;    	 	   									//Initial pressure head
    volScalarField hM = h;    	 	   									//Previous Picard Iteration level pressure head
    volScalarField hN = h;    	 	   									//Previous time level pressure head
    volScalarField hNm1 = h;  	 	   									//Previous to previous time level pressure head
    
	//------------------------------------------------------------------------------------------------------------------------------//
		
	#include "initializeSoilWaterParameters.H"							//Initializing Soil-water Characteristic Parameters
	
	//------------------------------------------------------------------------------------------------------------------------------//
	
	//Selection of the Convergence Criterion
	scalar Pic_tol = 0.0;
	
	if (convergenceCriterion.match("h-based"))							//Standard pressure head-based
	{
		Pic_tol = delta_htol;
	}
	else if (convergenceCriterion.match("mixed_h-based"))				//Mixed pressure head-based
	{
		Pic_tol = delta_htol; 
	}
	else if (convergenceCriterion.match("theta-based"))					//Standard moisture content-based
	{
		Pic_tol = delta_thetatol; 
	}
	
	//------------------------------------------------------------------------------------------------------------------------------//
	
	//Opening files for writing the MBE and Number of Picard Iterations for every time-step in a CSV file
	std::ofstream file1;
	file1.open("massBalanceError.csv");
	
	std::ofstream file2;
	file2.open("PicardIterations.csv");
	
	//------------------------------------------------------------------------------------------------------------------------------//
	
	//Time-Loop
    while (runTime.loop())
    {
		#include "readTimeControls.H"
		
		Info << "deltaT = " <<  runTime.deltaT().value() << endl;
        
        Info << "\nTime = " << runTime.timeName() << endl;	
        
        //--------------------------------------------------------------------------------------------------------------------------//
        
		//Picard Iteration Loop
		for(label count = 1; count < NPmax + 1; count++)
		{
			if (RicForm.match("h-based"))
			{
				#include "h-basedRE.H"
			}
			else if (RicForm.match("mixed"))
			{
				#include "mixedRE.H"
			}
			
			//----------------------------------------------------------------------------------------------------------------------//	
			
			#include "updateSoilWaterParameters.H" 						//Updating Soil-water Characteristic Parameters
			
			//----------------------------------------------------------------------------------------------------------------------//
								
			#include "convergenceCriteria.H" 							//Calculating the Convergence Criterion
				
			//----------------------------------------------------------------------------------------------------------------------//
			
			//Termination criterion for Picard Iteration Loop
			reduce(count, maxOp<label>());
			
			if ((error < Pic_tol) && (count >= NPmin))
			{
				NP = count;
				Info << "\nError = " << error << endl;
				Info << "\nNo.of Picard Iterations = " << NP << endl;
				break;
			}
			else if (count == NPmax)
			{
				NP = count;
				Info << "\nConvergence failed!" << endl;
				Info << "\nError = " << error << endl;
				Info << "\nNo.of Picard Iterations = " << NP << endl;
				
				runTime.setTime											//Set Runtime to previous time-level
				(
					oldTime, oldTimeIndex
				);
			 
				h = hN;													//Resetting 'h' and 'hN' field values
				hN = hNm1;
				
				thetaN = thetaNm1;										//Resetting 'theta' and 'thetaN' values
				
				#include "updateSoilWaterParameters.H" 					//Updating Soil-water Characteristic Parameters
				
				nMARK = 1;	
			}
			else
			{
				NP = count;
			}
	
			//----------------------------------------------------------------------------------------------------------------------//
			
			hM = h;
			thetaM = theta;
			
		} 																//End of Picard Iteration Loop
		
		//--------------------------------------------------------------------------------------------------------------------------//
		
		U = -K & (fvc::grad(h) + grad_z);								//Updating the cell-centred and face-centred flux
		phi = fvc::interpolate(U) & mesh.Sf();
		
		//--------------------------------------------------------------------------------------------------------------------------//
		
		#include "massBalanceError.H" 									//Calculation of Mass Balance Error (MBE)
		
		if (nMARK != 1)
		{			
			file1 << runTime.timeName() << "," << MB1 << "," << MB2 << "," << MBE << nl;	//Writing the MBE
			file2 << runTime.timeName() << "," << NP << nl;					//Writing the Number of Picard Iterations	
		}
		
		//--------------------------------------------------------------------------------------------------------------------------//
																		
		runTime.write();
		
		oldTime = runTime.value();
		oldTimeIndex = runTime.timeIndex();
		
		hNm1 = hN;
		hN = h;
		thetaNm1 = thetaN;
		thetaN = theta;
		nMARK = 0;
		
		#include "setDeltaT.H"											//Adjust time-step
		
	} 																	//End of Time-Loop
	
	//------------------------------------------------------------------------------------------------------------------------------//
	
	file1.close();
	file2.close();

	//------------------------------------------------------------------------------------------------------------------------------//
	
	std::clock_t endT= std::clock(); 									//End Time
	
	//Writing the Simulation Time
	scalar simTime;
	simTime = (endT - startT)/ (double) CLOCKS_PER_SEC;
	Info<< "\nCPU simulation time = " << simTime << nl << endl;
    
    Info<< "End\n" << endl;

    return 0;
}
