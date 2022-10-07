/*
tmp<volTensorField> tgradU(fvc::grad(U));
const volTensorField& gradU = tgradU();
*/
const word& wallName(runTime.controlDict().lookup("reactingWall"));
Info << "Reacting surface name: " << wallName << endl;
scalar loopNum = readScalar(runTime.controlDict().lookup("loopNum"));
scalar pseudoDt = readScalar(runTime.controlDict().lookup("pseudoDt"));
Info << "loopNum: " << loopNum << " pseudoDt : " << pseudoDt << nl << endl;

// ========== Diffusion ============ //

volScalarField sr
(
 "sr",
 Foam::sqrt(2.0) * mag(symm(fvc::grad(U)))
 );

volScalarField E_D_PLT
(
 "E_D_PLT",
D_PLT + 0.18*rRBC*rRBC*sr
 );

volScalarField E_D_ADP
(
 "E_D_ADP",
D_ADP + 0.18*rRBC*rRBC*sr
 );
volScalarField E_D_TXA2
(
 "E_D_TXA2",
D_TXA2 + 0.18*rRBC*rRBC*sr
 );
volScalarField E_D_PT
(
 "E_D_PT",
D_PT + 0.18*rRBC*rRBC*sr
 );
volScalarField E_D_TB
(
 "E_D_TB",
D_TB + 0.18*rRBC*rRBC*sr
 );
volScalarField E_D_ATI
(
 "E_D_AT",
D_ATI + 0.18*rRBC*rRBC*sr
 );

// Gamma
volScalarField Gamma
(
	"Gamma",
	k_1_T * H * ATI / (alpha * K_AT * K_T + alpha * K_AT * beta * TB + ATI * beta * TB)
);

volScalarField k_PA
(
    "k_PA",
    (w_adp * ADP / C_ADP + w_txa2 * TXA2 / C_TXA2 + w_tb * TB / C_TB) / t_act
);

scalar dT = runTime.deltaT().value();
runTime.setDeltaT(pseudoDt);
for (int loop = 0; loop < loopNum; loop++)
{

// ================================ //
//Info << "Calculation of Coeff D Done" << endl;

// Mesh
label wallId = mesh.boundaryMesh().findPatchID(wallName);
const vectorField& wallCenter = mesh.boundary()[wallId].Cf();
const vectorField& wallArea = mesh.boundary()[wallId].Sf();
const labelUList& faceCells = mesh.boundaryMesh()[wallId].faceCells();

// scalars
scalarField sr_Field = sr.boundaryField()[wallId];
scalarField& cThrom_Field = cellThrom.boundaryFieldRef()[wallId];
scalarField& fThrom_Field = faceThrom.boundaryFieldRef()[wallId];

scalarField D_PLT_Field = E_D_PLT.boundaryField()[wallId];
scalarField RPLT_Field = RPLT.boundaryField()[wallId];
scalarField APLT_Field = APLT.boundaryField()[wallId];

scalarField ADP_Field = ADP.boundaryField()[wallId];
scalarField TXA2_Field = TXA2.boundaryField()[wallId];
scalarField PT_Field = PT.boundaryField()[wallId];
scalarField TB_Field = TB.boundaryField()[wallId];

// reactions
scalarField& kpa_Field = k_PA.boundaryFieldRef()[wallId];

scalarField& m_Field = M.boundaryFieldRef()[wallId];
scalarField& mas_Field = M_as.boundaryFieldRef()[wallId];
scalarField& mr_Field = M_r.boundaryFieldRef()[wallId];
scalarField& mat_Field = M_at.boundaryFieldRef()[wallId];

scalarField& jr_Field = Jr.boundaryFieldRef()[wallId];
scalarField& ja_Field = Ja.boundaryFieldRef()[wallId];
scalarField& jpr_Field = Jpr.boundaryFieldRef()[wallId];
scalarField& jps_Field = Jps.boundaryFieldRef()[wallId];
scalarField& jpt_Field = Jpt.boundaryFieldRef()[wallId];
scalarField& jtb_Field = Jtb.boundaryFieldRef()[wallId];

// Info << "========Boundary==========" << endl;
std::queue<scalar> activeCellList;
std::set<scalar> cacheSet;
std::set<scalar>::iterator cacheIter;

std::queue<scalar> activeCCList;
std::set<scalar> cacheCCSet;
std::set<scalar>:: iterator ccIter;

forAll(wallCenter, faceI){
    // vector faceCenter = wallCenter[faceI];
    // Pre Values
    dimensionedScalar k_pa = kpa_Field[faceI] * kUnit * (kpa_Field[faceI] >= 1);
    
    // values are given
    /* 
    dimensionedScalar k_as = Sh * D_PLT_Field[faceI] / cL * dUnit;
    dimensionedScalar k_rs = Sh * D_PLT_Field[faceI] / cL * dUnit;
    dimensionedScalar k_aa = Sh * D_PLT_Field[faceI] / cL * dUnit;
    */
    // surface coverage between 0 and 1
    // dimensionedScalar s = 1. - m_Field[faceI] / M_inf * mUnit; // Eq 11
    
    scalar s = fmax(fmin(1. - m_Field[faceI] / M_inf.value(), 1.), 0);
    if (faceI % 1000 == 0)
        Info << "SS: " <<  s << endl;
    // surface reaction unit: Quantity / m2 / t. temporal flux
    dimensionedScalar jr = s * k_rs * max(RPLT_Field[faceI], 0) * cUnit;
    dimensionedScalar ja = (s * k_as + mas_Field[faceI] * mUnit / M_inf * k_aa) * APLT_Field[faceI] * cUnit;
    dimensionedScalar jpr = -lmd_ADP * (theta * s * k_rs * RPLT_Field[faceI] * cUnit + k_pa * mr_Field[faceI] * mUnit);
    dimensionedScalar jps = -mat_Field[faceI] * s_TXA2 * mUnit;
    dimensionedScalar jpt = beta * (mat_Field[faceI] * phi_APLT + mr_Field[faceI] * phi_RPLT) * PT_Field[faceI] * mUnit * cUnit;
    dimensionedScalar jtb = -(mat_Field[faceI] * phi_APLT + mr_Field[faceI] * phi_RPLT) * PT_Field[faceI] * mUnit * cUnit;

    dimensionedScalar dm = jr + k_as * s * APLT_Field[faceI] * cUnit;
    //dimensionedScalar dm = jr_Field[faceI] * jUnit + k_as * s * APLT_Field[faceI] * cUnit;
    dimensionedScalar m = m_Field[faceI] *mUnit + dm * runTime.deltaT(); // Eq 12
    
    dimensionedScalar dmas = theta * jr + k_as * s * APLT_Field[faceI] * cUnit + k_pa * mr_Field[faceI] * mUnit;
    // dimensionedScalar dmas = theta * jr_Field[faceI] * jUnit + k_as * s * APLT_Field[faceI] * cUnit + k_pa * mr_Field[faceI] * mUnit;
    dimensionedScalar mas = mas_Field[faceI] * mUnit + dmas * runTime.deltaT(); // Eq 13
    
    dimensionedScalar dmr = (1. - theta) * jr - k_pa * mr_Field[faceI] * mUnit;
    // dimensionedScalar dmr = (1. - theta) * jr_Field[faceI] * jUnit - k_pa * mr_Field[faceI] * mUnit;
    dimensionedScalar mr = mr_Field[faceI] * mUnit + dmr * runTime.deltaT(); // Eq 14
    
    dimensionedScalar dmat = theta * jr + ja + k_pa * mr_Field[faceI] * mUnit;
    // dimensionedScalar dmat = theta * jr_Field[faceI] * jUnit + ja_Field[faceI] * jUnit + k_pa * mr_Field[faceI] * mUnit;
    dimensionedScalar mat = mat_Field[faceI] * mUnit + dmat * runTime.deltaT(); // Eq 15
    
    // New Values
    m_Field[faceI] = m.value();
    mas_Field[faceI] = mas.value();
    mr_Field[faceI] = mr.value();
    mat_Field[faceI] = mat.value();

    jr_Field[faceI] = jr.value();// * runTime.deltaT().value();
    ja_Field[faceI] = max(0, ja.value());// * runTime.deltaT().value();
    jpr_Field[faceI] = jpr.value();// * runTime.deltaT().value();
    jps_Field[faceI] = jps.value();// * runTime.deltaT().value();
    jpt_Field[faceI] = jpt.value();// * runTime.deltaT().value();
    jtb_Field[faceI] = jtb.value();// * runTime.deltaT().value();

    // Platelet embolization model
    //scalar F_emb = exp(-0.0095 * E_D_PLT[faceCells[faceI]] * sr_Field[faceI]);

    // Platelet plug growth model
    scalar annex = (jr_Field[faceI] + ja_Field[faceI]) * (1 + 0.48) * v_plt.value() * runTime.deltaT().value() * ampPlt.value(); // mPLT
    cThrom_Field[faceI] += annex; // Accumulate distance value
    label faceOwner = faceCells[faceI];
    scalar dis = mag(mesh.C()[faceOwner] - wallCenter[faceI]);
    // cell find the distance from the wall

    if (cThrom_Field[faceI] - dis > 0)
    {
        cellThrom[faceOwner] = cThrom_Field[faceI] - dis;
        activeCellList.push(faceOwner);
    }
    
    // Annex cells based on volume
    // density of cells per cm^3
    scalar generation = annex * mag(wallArea[faceI]);
    Info << faceI << " " << faceOwner << " " << generation << endl;
    Info << "face area: " << mag(wallArea[faceI]) << " cell v: " << mesh.V()[faceOwner] << endl;
    Info << "sum: " << (cThrom_Field[faceI] * mag(wallArea[faceI])) << endl;
    faceThrom[faceOwner] += generation;
    fThrom_Field[faceI] += generation;
    if (faceThrom[faceOwner] > mesh.V()[faceOwner])
    {
        // faceThrom[faceOwner] = generation - mesh.V()[faceOwner];
        activeCCList.push(faceOwner);
    }
    Info << "faceThrom: " << faceThrom[faceOwner] << endl;
}

Info << "Face Owner: " << activeCellList.size() << endl;
while (activeCellList.size() > 0){
    int num = activeCellList.size();
    cacheSet.clear();
    for(int ci_ = 0; ci_ < num; ci_++){
        scalar centerCellId = activeCellList.front();
        activeCellList.pop();
        scalar centerR = cellThrom[centerCellId];
        vector centerX = mesh.C()[centerCellId];

        const labelUList& neighbours = mesh.cellCells()[centerCellId];
        forAll(neighbours, cc_){
            scalar curCellId = neighbours[cc_];
            vector cX = mesh.C()[curCellId];
            scalar cR = mag(cX - centerX);
            // Info << "id: " << curCellId << " T: " << Throm[curCellId] << " X: " << cX << " dis: " << cR << endl;
            // Info << "centerR: " << centerR << " cR: " << cR << " Throm: " << Throm[curCellId] << endl;

            if (centerR - cR > cellThrom[curCellId])
            {
                cellThrom[curCellId] = centerR - cR;
                if (cacheSet.find(curCellId) == cacheSet.end())
                {
                    activeCellList.push(curCellId);
                    cacheSet.insert(curCellId);
                }
            }
            /*
            if ((cacheSet.find(curCellId) == cacheSet.end()) && Throm[curCellId] > 0)
            {
                activeCellList.push(curCellId);
                cacheSet.insert(curCellId);
            }
            */
        }
    }
    // Info << "Next Round: " << activeCellList.size() << endl;
}

const labelUList& owner = mesh.owner();
const labelUList& neighbour = mesh.neighbour();

Info << "Face Owner: " << activeCCList.size() << endl;
while (activeCCList.size()){
    int num = activeCCList.size();
    cacheCCSet.clear();
    for (int ci_ = 0; ci_ < num; ci_++){
        scalar centerCellId = activeCCList.front();
        activeCCList.pop();
        scalar centerQ = faceThrom[centerCellId];
        scalar centerV = mesh.V()[centerCellId];
        scalar centerDensity = centerQ / centerV;
        Info << "Cell Id: " << centerCellId << endl;

        const labelUList& neighbourFaces = mesh.cells()[centerCellId];
        scalar totalArea = 0;
        std::list<std::pair<scalar, scalar>> activeFaces;
        forAll(neighbourFaces, ff_){
            scalar faceId_ = neighbourFaces[ff_];
            scalar left_ = mesh.owner()[faceId_];
            scalar right_ = mesh.neighbour()[faceId_];
            //Info << "left: " << left_ << " " << right_ << " " << centerCellId << endl;
            scalar area_ = mesh.magSf()[faceId_];
            scalar cellId_ = (faceId_ == left_)? left_ : right_;
            scalar thromV_ = faceThrom[cellId_];
            scalar meshV_ = mesh.V()[cellId_];
            scalar density_ = thromV_ / meshV_;

            if ((density_ < centerDensity) && (meshV_ > 0))
            {
                totalArea += area_;
                activeFaces.push_back({cellId_, area_});
            }
        
        }

        for(auto ff_ = activeFaces.begin(); ff_ != activeFaces.end();ff_++){
            scalar cellId_ = ff_->first;
            scalar frac_ = ff_->second / totalArea;
            faceThrom[cellId_] += (faceThrom[centerCellId] - mesh.V()[centerCellId]) * frac_;
            if ((faceThrom[cellId_] > mesh.V()[cellId_]) && (cacheCCSet.find(cellId_) == cacheCCSet.end()))
            {
                activeCCList.push(cellId_);
                cacheCCSet.insert(cellId_);
            }
        }
    }
}


// Info << "==========Build Matrix=========" << endl;
// SPECIES 
//Info << "Build RPLT" << endl;
fvScalarMatrix RPLTEqn // === unit: PLT / ml === //
(
    fvm::ddt(RPLT) + fvm::div(phi, RPLT) - fvm::laplacian(E_D_PLT, RPLT)
    == -k_PA * RPLT // source term
 );
//Info << "Build APLT" << endl;
 fvScalarMatrix APLTEqn // === unit: PLT / ml === //
 (
     fvm::ddt(APLT) + fvm::div(phi, APLT) - fvm::laplacian(E_D_PLT, APLT)
     == k_PA * RPLT// source term
 );
//Info << "Build ADP" << endl;
fvScalarMatrix ADPEqn  // === unit: muM (nmol / ml) === //
(
     fvm::ddt(ADP) + fvm::div(phi, ADP) - fvm::laplacian(E_D_ADP, ADP)
     == lmd_ADP * k_PA * RPLT 
);
//Info << "Build TXA2" << endl;
fvScalarMatrix TXA2Eqn
(
    fvm::ddt(TXA2) + fvm::div(phi, TXA2) - fvm::laplacian(E_D_TXA2, TXA2)
    == s_TXA2 * APLT - fvm::Sp(k_TXA2, TXA2) //k_TXA2 * TXA2
);
//Info << "Build PT" << endl;
fvScalarMatrix PTEqn  // === unit: muM (nmol / ml) === //
(
    fvm::ddt(PT) + fvm::div(phi, PT) - fvm::laplacian(E_D_PT, PT)
    == -fvm::Sp(beta * (phi_APLT * APLT + phi_RPLT * RPLT), PT)//-beta * (phi_APLT * APLT + phi_RPLT * RPLT) * PT
);
//Info << "==========Calculate Matrix==========" << endl;
APLTEqn.solve();
RPLTEqn.solve();
ADPEqn.solve();
TXA2Eqn.solve();
PTEqn.solve();

if (H.value() > 0)
    {
        Info << "Build TB with H" << endl;
        fvScalarMatrix TBEqn // === unit: U / ml === //
        (
            fvm::ddt(TB) + fvm::div(phi, TB) - fvm::laplacian(E_D_TB, TB)
            //== -Gamma * TB + PT * (phi_APLT * APLT + phi_RPLT * RPLT)
            == -fvm::Sp(Gamma, TB) + PT * (phi_APLT * APLT + phi_RPLT * RPLT)
        );
        //Info << "Build ATI" << endl;
        fvScalarMatrix ATIEqn //=== unit: muM (nmol / ml) === //
        (
            fvm::ddt(ATI) + fvm::div(phi, ATI) - fvm::laplacian(E_D_ATI, ATI)
            == -Gamma * beta * TB
        );

        TBEqn.solve();
        ATIEqn.solve();
    }
else
    {
        Info << "Build TB without H" << endl;
        fvScalarMatrix TBEqn // === unit: U / ml === //
        (
            fvm::ddt(TB) + fvm::div(phi, TB) - fvm::laplacian(E_D_TB, TB)
            //== -k2 * ATI * TB + PT * (phi_APLT * APLT + phi_RPLT * RPLT)
            == -fvm::Sp(k2 * ATI, TB) + PT * (phi_APLT * APLT + phi_RPLT * RPLT)
        );
        //Info << "Build ATI" << endl;
        fvScalarMatrix ATIEqn //=== unit: muM (nmol / ml) === //
        (
            fvm::ddt(ATI) + fvm::div(phi, ATI) - fvm::laplacian(E_D_ATI, ATI)
            // == -k2 * ATI * beta * TB
            == -fvm::Sp(k2 * beta * TB, ATI)
        );

        TBEqn.solve();
        ATIEqn.solve();
    }

}

// Calculate gradient of APLT

//volVectorField apltFlux ("apltFlux", fvc::grad(APLT));

// Update Velocity source term
volScalarField nu = laminarTransport.nu();
dimensionedScalar k_plug = Foam::pow(0.48, 3) / (36. * 5 * Foam::pow(1 - 0.48, 2)) * Foam::pow(r_plt, 2);

forAll(mesh.C(), cI){
    sU[cI] = (cellThrom[cI] > 0) * nu[cI] / k_plug.value();
}



runTime.setDeltaT(dT);

if (runTime.outputTime())
    {
        k_PA.write();
        //Gamma.write();
        sr.write();
        //E_D_PLT.write();
        //apltFlux.write();
    }
