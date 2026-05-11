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

Description
    Calculates species mass fractions and thermodynamic properties
    from given Z, varZ , PV and varPV fields

    @author Likun Ma, Delft University of Technology
    @email  malikun-2005@hotmail.com
    @version 14.11.2014
    @version 30.01.2020 mb

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "tableSolver.H"
#include "PVtableSolver.H"

#include "fvCFD.H"
#include "CombustionModel.H"  //changed by senbin
#include "psiReactionThermo.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be reconstructed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    //reading thermo properties, added by senbin
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    const IOdictionary combProps
    (
        IOobject
        (
            "combustionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const IOdictionary tableProps
    (
        IOobject
        (
            "tableProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const IOdictionary PVtableProps
    (
        IOobject
        (
            "PVtableProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const word combTypeName = combProps.lookup("combustionModel");
    const label tempOpen = combTypeName.find('<');
    const word modelType = combTypeName(0, tempOpen);
    dictionary coeffs_(combProps.subDict(modelType + "Coeffs"));

    Switch useProgressVariableVariance_(coeffs_.lookup("useProgressVariableVariance"));

    //- 2D table for minimum and maximum progress variable
    scalarList Params2d(2, 0.0);
    
    //- 4D independent variables (Z, varZ, PV, varPV)
    scalarList Params4d(4, 0.0);
    
    double Zvar, Cvar;


    List<List<int> > ubParams4dIF(mesh.cells().size()), ubParams2dIF(mesh.cells().size());
    List<scalarList> posParams4dIF(mesh.cells().size()), posParams2dIF(mesh.cells().size());
    List<List<int> > ubParams4dBF(mesh.faces().size()), ubParams2dBF(mesh.faces().size());
    List<scalarList> posParams4dBF(mesh.faces().size()), posParams2dBF(mesh.faces().size());

    //- Minimum and Maximum unscaled progress variable
    scalar YcMinCells, YcMaxCells, YcMinPatchs, YcMaxPatchs;

    //- Scaled progress variable   
    scalar Cmean;
    scalar Cfluc(0.0);

    //- Scale progress variable and its variance
    scalar fc, gc, hc;
    scalar Yu2I, YuYbI, Yb2I;

    wordList tableNames(thermo.composition().species());
    tableNames.append("T");   //--- Read 'T' table
    tableNames.append("SourceYc");
    Foam::combustionModels::tableSolver solver(Foam::combustionModels::tableSolver(mesh, tableNames));
    
    hashedWordList PVtableNames;
    PVtableNames.clear();
    PVtableNames.append("Ycmin");
    PVtableNames.append("Ycmax");
    if (useProgressVariableVariance_)
    {
        PVtableNames.append("Yu2I");
        PVtableNames.append("YuYbI");
        PVtableNames.append("Yb2I");
    }

    Foam::combustionModels::PVtableSolver PVsolver(Foam::combustionModels::PVtableSolver(mesh, PVtableNames));
    
    
    PtrList<volScalarField>& Y(thermo.composition().Y());
//    volScalarField& hs(thermo.he());   // Changed l.Ma, 25-06-2014    

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << nl << "Time = " << runTime.timeName() << nl << endl;

        volScalarField ZmeanCells
        (
            IOobject
            (
                "Zmean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField ZflucCells
        (
            IOobject
            (
                "Zfluc",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField YcmeanCells
        (
            IOobject
            (
                "Ycmean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField YcflucCells
        (
            IOobject
            (
                "Ycfluc",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

       // Interpolate for internal Field
       forAll(Y, i)
       {
            scalarField& YCells = Y[i].ref();   //internalField(); senbin

            forAll(ZmeanCells, cellI)
            {
                if (i == 0)  // Zeta and scaledPV determined once is enough
                {
            
                    //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as an independent variable
                    Zvar = ZflucCells[cellI]/max(ZmeanCells[cellI]*(1 - ZmeanCells[cellI]), SMALL);

                    //- Calculate scaled progress variable   
                    Params2d[0] = max(min(Zvar, 0.99), 0.0);
                    Params2d[1] = max(min(ZmeanCells[cellI], 1.0), 0.0);
            
                    ubParams2dIF[cellI] = PVsolver.upperBounds(Params2d);
                    posParams2dIF[cellI] = PVsolver.position(ubParams2dIF[cellI], Params2d);
        
                    //- Update minimum and maximum unscaled progress variable
                    YcMinCells = PVsolver.interpolate(ubParams2dIF[cellI], posParams2dIF[cellI], 0);
                    YcMaxCells = PVsolver.interpolate(ubParams2dIF[cellI], posParams2dIF[cellI], 1);
        
                    //- Update scalaed progress variable
                    Cmean = (YcmeanCells[cellI] - YcMinCells)/max((YcMaxCells - YcMinCells), SMALL);

                    //- How about the progress variable variance?
                    if (useProgressVariableVariance_)
                    {
                        //- calculate scaled progress variable variance from unscaled propress variance
                        //- Reference 1: A progress variable approach based on premixed flamelets for turbulent combustion modeling, B.A. Albrecht, W.J.S.Ramaekers et al
                        //- Eq.(10) and Eq.(17)
                        //- Reference 2: A premixed Flamelet-PDF Model for Biomass Combustion in a Grate Furnace, Energy & Fuels, 2008, Albrecht, Oijen et al
                        //- Eq.(10) and Eq.(11)
                        
                        Yu2I = PVsolver.interpolate(ubParams2dIF[cellI], posParams2dIF[cellI], 2);
                        YuYbI = PVsolver.interpolate(ubParams2dIF[cellI], posParams2dIF[cellI], 3);
                        Yb2I = PVsolver.interpolate(ubParams2dIF[cellI], posParams2dIF[cellI], 4);
                        
                        fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2*YuYbI + Yu2I;
                        Cfluc = (YcflucCells[cellI] + sqr(YcmeanCells[cellI]) - fc - 2*gc*Cmean)/max(hc, SMALL) - sqr(Cmean); 
                    }
                    else
                    {
                        Cfluc = 0.0;
                    }

                    Cvar = Cfluc/max(Cmean*(1 - Cmean), SMALL);

                    Params4d[0] = max(min(Cvar, 0.99), 0.0);
                    Params4d[1] = max(min(Cmean, 1.0), 0.0);
                    Params4d[2] = max(min(Zvar, 0.99), 0.0);
                    Params4d[3] = max(min(ZmeanCells[cellI], 1.0), 0.0);

                    //- find up-bound and pos for table interpolation
                    ubParams4dIF[cellI] = solver.upperBounds(Params4d);
                    posParams4dIF[cellI] = solver.position(ubParams4dIF[cellI], Params4d);

                    // YCells[cellI] = solver.interpolate(ubParams4dIF[cellI], posParams4dIF[cellI], i);
                }

                //- Update species
                YCells[cellI] = solver.interpolate(ubParams4dIF[cellI], posParams4dIF[cellI], i); 
            }
       }

       // Interpolate for patches
       forAll(ZmeanCells.boundaryField(), patchi)    // Changed L.Ma, 24-06-2014
       {
            const fvPatchScalarField& pYcfluc = YcflucCells.boundaryField()[patchi];
            const fvPatchScalarField& pYcmean = YcmeanCells.boundaryField()[patchi];
            const fvPatchScalarField& pZfluc = ZflucCells.boundaryField()[patchi];
            const fvPatchScalarField& pZmean = ZmeanCells.boundaryField()[patchi];

            forAll(Y, i)
            {
                fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchi];

                forAll(pY , facei)
                {
                    if (i == 0)
                    {
                        Zvar = pZfluc[facei]/max(pZmean[facei]*(1 - pZmean[facei]), SMALL);
                        //- Calculate scaled progress variable 
                        Params2d[0] = max(min(Zvar, 0.99), 0.0);
                        Params2d[1] = max(min(pZmean[facei], 1.0), 0.0);
        
                        ubParams2dBF[facei] = PVsolver.upperBounds(Params2d);
                        posParams2dBF[facei] = PVsolver.position(ubParams2dBF[facei], Params2d);
        
                        //- Update minimum and maximum unscaled progress variable
                        YcMinPatchs = PVsolver.interpolate(ubParams2dBF[facei], posParams2dBF[facei], 0);
                        YcMaxPatchs = PVsolver.interpolate(ubParams2dBF[facei], posParams2dBF[facei], 1);
        
                        //- Update scalaed progress variable
                        Cmean = (pYcmean[facei] - YcMinPatchs)/max((YcMaxPatchs - YcMinPatchs), SMALL);

                        //- How about the progress variable variance?
                        if (useProgressVariableVariance_)
                        {    
                            Yu2I = PVsolver.interpolate(ubParams2dBF[facei], posParams2dBF[facei], 2);
                            YuYbI = PVsolver.interpolate(ubParams2dBF[facei], posParams2dBF[facei], 3);
                            Yb2I = PVsolver.interpolate(ubParams2dBF[facei], posParams2dBF[facei], 4);
                            
                            fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2*YuYbI + Yu2I;
                            Cfluc =  (pYcfluc[facei] + sqr(pYcmean[facei]) - fc - 2*gc*Cmean)/max(hc, SMALL) - sqr(Cmean); 
                        }
                        else
                        {
                            Cfluc = 0.0;
                        }

                        Cvar = Cfluc/max(Cmean*(1 - Cmean), SMALL);

                        Params4d[0] = max(min(Cvar, 1.0), 0.0);
                        Params4d[1] = max(min(Cmean, 1.0), 0.0);
                        Params4d[2] = max(min(Zvar, 0.99), 0.0);
                        Params4d[3] = max(min(pZmean[facei], 1.0), 0.0);
                
                        ubParams4dBF[facei] = solver.upperBounds(Params4d);
                        posParams4dBF[facei] = solver.position(ubParams4dBF[facei], Params4d);

                        // pY[facei] = solver.interpolate(ubParams4dBF[facei], posParams4dBF[facei], i);
                    }
                    //- update speces
                    pY[facei] = solver.interpolate(ubParams4dBF[facei], posParams4dBF[facei], i);
                }
            }
       }
  

        if (selectedFields.empty())
        {
        	forAll(Y, i)
            {
               Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        	   thermo.composition().Y()[i].write();
            }
        }
        else
        {
        	forAll(Y, i)
            {
        	   if (selectedFields[thermo.composition().Y()[i].name()])
        	   {
                   Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        		   thermo.composition().Y()[i].write();
        	   }
            }
        }

    }

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "End\n" << endl;

	return 0;
}
