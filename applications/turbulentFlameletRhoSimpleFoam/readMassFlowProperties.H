// Initialisierung am Anfang vor dem Loop!, danach werden diese Zeilen nicht mehr durchlaufen!
// Beschreibung im write Mass Flow

Info<< "Reading mass flow properties\n" << endl;

IOdictionary massFlowDict
(
    IOobject
    (
        "massFlowProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

wordList inletPatches(massFlowDict.lookup("inletPatches"));	// Hole die Oberflächen für Inlet - Liste wird vom User vorgegeben
wordList outletPatches(massFlowDict.lookup("outletPatches"));	// Äquivalent zu oben

fileName path = runTime.caseName();				

ofstream outInfo("outInfo.case", ios::out);			// Name des Files <outInfo.Case>

outInfo << "time   "  						// Erste Zeile im File	!
        << "m_in   "  
        << "v_in   "  
        << "H_in   "  
        << "m_out  " 
        << "v_out  " 
        << "H_out  " 
        << "m_net  " 
        << "v_net  " 
        << "H_net  " 
        << endl;

