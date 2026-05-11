#include "createRDeltaT.H"

Info<< "Creating combustion model\n" << endl;

// The initialization of thermo changed in of6
// Thermo model
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();

const volScalarField& psi = thermo.psi();

// Pressure
volScalarField& p = thermo.p();

// Temperature
volScalarField& T = thermo.T();

// Density
// volScalarField rho
// (
//     IOobject
//     (
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     thermo.rho()
// );
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// Velocity
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);
mesh.setFluxRequired(p.name());


Info << "Creating turbulence model.\n" << nl;

// Set some booleans for the turbulence models 
bool laminarBool  = true;
bool RASBool      = false;
// bool kEpsilonBool = false;  // L
bool LESBool      = false;

// List of all available RAS models  // L
wordList RASModels = {"kEpsilon",
                      "kOmega",
                      "kOmegaSST",
                      "kOmegaSSTLM",
                      "kOmegaSSTSAS",
                      "LaunderSharmaKE",
                      "LRR",
                      "realizableKE",
                      "RNGkEpsilon",
                      "SpalartAllmaras",
                      "SSG",
                      "v2f"};

// List of all available LES models
wordList LESModels = {"DeardorffDiffStress",
                      "Smagorinsky",
                      "SpalartAllmarasDDES",
                      "SpalartAllmarasDES",
                      "SpalartAllmarasIDDES",
                      "WALE",
                      "dynamicKEqn",
                      "dynamicLagrangian",
                      "kEqn",
                      "kOmegaSSTDDES",
                      "kOmegaSSTDES",
                      "kOmegaSSTIDDES"};

autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

const word turbModelName(turbulence->type());
Info<< "Turbulent Model Used: " << turbModelName << endl;

if (turbulence->type() != "Stokes"){ // If not laminar
    laminarBool = false;
    // if (turbulence->type() == "kEpsilon"){ // If RANS and kEpsilon
        // kEpsilonBool = true;
    if (std::find(std::begin(RASModels), std::end(RASModels), turbulence->type()) != std::end(RASModels)){
        RASBool = true;
    }
    else if (std::find(std::begin(LESModels), std::end(LESModels), turbulence->type()) != std::end(LESModels)){ // If LES
        LESBool = true;
    }
    else { // Only the RANS kEpsilon and all LES models work with this solver
        FatalErrorInFunction << "\nPlease provide valid turbulence model. Available models are:\n\nRANS:\tkEpsilon\nLES:\tall models" << exit(FatalError);
    }
}

// Mixture fraction
volScalarField& Zmean = turbulence().Zmean();
// Mixture fraction var
volScalarField& Zfluc = turbulence().Zfluc();

// Progress variable
volScalarField& Ycmean = turbulence().Ycmean();
// Progress variable var
volScalarField& Ycfluc = turbulence().Ycfluc();

// Create interpolation scheme for control variables
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

fields.add(Zmean);
fields.add(Zfluc);
fields.add(Ycmean);
fields.add(Ycfluc);


Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);
const word combModelName(reaction->type());
Info<< "Combustion Model Used: " << combModelName << endl;

// Set the turbulence into the combustion model, exist in of2 and of4, but not exist in of6
//reaction->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvOptions.H"


// Turbulent scalar dissipation of the progress variable
volScalarField XiPV
(
    IOobject
    (
        "XiPV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,0,-1,0,0,0,0), Zero)
);

volScalarField XiZ
(
    IOobject
    (
        "XiZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,0,-1,0,0,0,0), Zero)
);

// Turbulent Schmidt number
// volScalarField Sct
// (
//     IOobject
//     (
//         "Sct",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar(dimless, 0.7)
// );

// LES turbulent scalar dissipation model parameter
volScalarField CEpsilon
(
    IOobject
    (
        "CEpsilon",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1.8)
);

// Thermal diffusivity
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.alpha()
);


// Arbritrary small k value to prevent division by zero
volScalarField smallk
(
    IOobject
    (
        "smallk",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,2,-2,0,0,0,0), 1E-10)
);

// Field with cell volumes. This is needed to ensure type coherence in the Xi computation
volScalarField cellVolumes
(
    IOobject
    (
        "cellVolumes",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,3,0,0,0,0,0), 1E-10)
);
