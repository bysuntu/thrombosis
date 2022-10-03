Info << "Reading transportProperties" << nl << endl;
// (g cm s 0 Concentration 0 0)
IOdictionary transportProperties
(
 IOobject
 (
  "transportProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

// Diffusivity Table 1 of Paper 1
// platelet
dimensionedScalar D_PLT(transportProperties.lookup("D_PLT")); 
// ADP
dimensionedScalar D_ADP(transportProperties.lookup("D_ADP"));
dimensionedScalar D_TXA2(transportProperties.lookup("D_TXA2"));
// prothombin
dimensionedScalar D_PT(transportProperties.lookup("D_PT"));
// thrombin
dimensionedScalar D_TB(transportProperties.lookup("D_TB"));
// atiii
dimensionedScalar D_ATI(transportProperties.lookup("D_ATI"));

// k
dimensionedScalar k_TXA2(transportProperties.lookup("k_TXA2"));

// lmb
dimensionedScalar lmd_ADP(transportProperties.lookup("lmd_ADP")); // nmol / PLT

// s
dimensionedScalar s_TXA2(transportProperties.lookup("s_TXA2")); // per 1e9 PLT

// a
dimensionedScalar C_ADP(transportProperties.lookup("C_ADP"));
dimensionedScalar C_TXA2(transportProperties.lookup("C_TXA2"));
dimensionedScalar C_TB(transportProperties.lookup("C_TB"));

dimensionedScalar w_adp(transportProperties.lookup("w_adp"));
dimensionedScalar w_txa2(transportProperties.lookup("w_txa2"));
dimensionedScalar w_tb(transportProperties.lookup("w_tb"));

dimensionedScalar k_1_T(transportProperties.lookup("k_1_T"));
dimensionedScalar K_AT(transportProperties.lookup("K_AT"));
dimensionedScalar K_T(transportProperties.lookup("K_T"));
dimensionedScalar alpha(transportProperties.lookup("alpha"));
dimensionedScalar beta(transportProperties.lookup("beta"));
dimensionedScalar M_inf(transportProperties.lookup("M_inf"));
dimensionedScalar phi_APLT(transportProperties.lookup("phi_APLT"));
dimensionedScalar phi_RPLT(transportProperties.lookup("phi_RPLT"));

// Heparin
dimensionedScalar H(transportProperties.lookup("H")); // Heparin
dimensionedScalar Sh(transportProperties.lookup("Sh")); // Sherwood Number
dimensionedScalar cL(transportProperties.lookup("cL")); // characteristic Length

dimensionedScalar theta(transportProperties.lookup("theta")); // Theta
dimensionedScalar k2(transportProperties.lookup("k2")); // k2

dimensionedScalar rRBC(transportProperties.lookup("rRBC")); // radius of RBC in cm

dimensionedScalar k_as(transportProperties.lookup("k_as")); // reaction rate cm / s. activate platlets and the surface.
dimensionedScalar k_rs(transportProperties.lookup("k_rs")); // reaction rate cm / s. resting platelets and the surface.
dimensionedScalar k_aa(transportProperties.lookup("k_aa")); // reaction rate cm / s. depositng and surface-bound activated platelets.

dimensionedScalar t_act("t_act", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.);

volScalarField RPLT
(
    IOobject
    (
        "RPLT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField APLT
(
    IOobject
    (
        "APLT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ADP
(
    IOobject
    (
        "ADP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField TXA2
(
    IOobject
    (
        "TXA2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField PT
(
    IOobject
    (
        "PT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField TB
(
    IOobject
    (
        "TB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField ATI
(
    IOobject
    (
        "ATI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField M 
{
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("M", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField M_as
{
    IOobject
    (
        "M_as",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("M_as", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField M_r
{
    IOobject
    (
        "M_r",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("M_r", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField M_at
{
    IOobject
    (
        "M_at",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("M_at", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField faceThrom
{
    IOobject
    (
        "faceThrom",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("faceThrom", dimensionSet(0, 0, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField cellThrom
{
    IOobject
    (
        "cellThrom",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("cellThrom", dimensionSet(0, 0, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
};

volScalarField Jr
(
    IOobject
    (
        "Jr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Jr", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
);

volScalarField Ja
(
    IOobject
    (
        "Ja",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ja", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
);

volScalarField Jpr
(
    IOobject
    (
        "Jpr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("Jpr", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
);

volScalarField Jps
(
    IOobject
    (
        "Jps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("Jps", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
);

volScalarField Jpt
(
    IOobject
    (
        "Jpt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("Jpt", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)

);

volScalarField Jtb
(
    IOobject
    (
        "Jtb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("Jtb", dimensionSet(0, -2, 0, 0, 1, 0, 0), pTraits<scalar>::zero)
);

volScalarField sU
(
    IOobject
    (
        "sU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


dimensionedScalar cUnit("cUnit", ADP.dimensions(), 1.);
dimensionedScalar dUnit("dUnit", D_PLT.dimensions(), 1.);
dimensionedScalar mUnit("mUnit", M.dimensions(), 1.);
dimensionedScalar jUnit("jUnit", Jr.dimensions(), 1.);
dimensionedScalar kUnit("kUnit", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.);


#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
