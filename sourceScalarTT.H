volScalarField sAPLT
(
    IOobject
    (
        "sAPLT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sAPLT", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sRPLT
(
    IOobject
    (
        "sRPLT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sRPLT", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sADP
(
    IOobject
    (
        "sADP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sADP", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sTXA2
(
    IOobject
    (
        "sTXA2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sTXA2", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sPT
(
    IOobject
    (
        "sPT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sPT", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sTB
(
    IOobject
    (
        "sTB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sTB", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);

volScalarField sATI
(
    IOobject
    (
        "sATI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sATI", dimensionSet(0, -3, -1, 0, 1, 0, 0), pTraits<scale>::zero)
);
