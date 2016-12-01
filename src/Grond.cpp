#include <Grond.h>


Grond::Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,float _drogeMassaDichtheid, std::string _Naam,double c_v,double k_s=0.0):
    samendrukkingsCoeff(_samendrukkingscoeff),bovengrens(_bovengrens), Naam(_Naam),
    ondergrens(_ondergrens),drogeMassDichtheid(_drogeMassaDichtheid),c_v(c_v){

    laagdikte = std::abs(bovengrens - ondergrens);
    gen_js();
    gen_msg();
}

Grond::Grond(json grondJS){

    samendrukkingsCoeff = grondJS["samendrukkingsCoeff"].get<double>();
    ondergrens = grondJS["ondergrens"].get<double>();
    bovengrens = grondJS["bovengrens"].get<double>();
    drogeMassDichtheid = grondJS["drogeMassaDichtheid"].get<double>();
    Naam = grondJS["Naam"].get<std::string>();
    c_v = (grondJS["C_v"].get<double>()<0.00000001)?0: grondJS["C_v"].get<double>();
    k_s= (grondJS["k_s"].get<double>()<0.00000001) ? 0 : grondJS["k_s"].get<double>();
    laagdikte = std::abs(bovengrens - ondergrens);
    ground_js = grondJS;
    //gen_js();
    gen_msg();
}

void Grond::gen_msg() {
    std::ostringstream str;

    str << "Grondlaag " << Naam << " met samendrukkingscoefficient C= " <<
        std::setprecision(decimalPrecisionInShout) << samendrukkingsCoeff <<
        " \n De laag gaat van " << std::setprecision(decimalPrecisionInShout) << (bovengrens)
        << "m tot " << std::setprecision(decimalPrecisionInShout) << (ondergrens) <<
        "m \n De massadichtheid van de grond bedraagt: " << std::setprecision(decimalPrecisionInShout)
        << (drogeMassDichtheid) << "kN/m^3\n" << "met C_v: " << c_v <<
        "m^2/s\n" << "met doorlatendheid: " << k_s <<"m/s\n";
    Message = str.str();
    str.clear();
}

double Grond::getLaagDikteNaPrim() {
    dikteNaPrim = laagdikte - primZetting;
    return dikteNaPrim;
}


void Grond::gen_js() {
    json js;
    js["samendrukkingsCoeff"] = samendrukkingsCoeff;
    js["ondergrens"] = ondergrens;
    js["bovengrens"] = bovengrens;
    js["drogeMassaDichtheid"] = drogeMassDichtheid;
    js["Naam"] = Naam;
    js["C_v"] = c_v;
    js["k_s"] = k_s;
    ground_js = js;
}

Grond::~Grond()
{
}

std::string Grond::shout()
{
    return Message;
}

Zettingsberekening::~Zettingsberekening()
{
}

Zettingsberekening::Zettingsberekening()
{
}

Zettingsberekening::Zettingsberekening(json js){

    belastingsType = BelastingsType(js["belasting"]);

    std::vector<json> grondlagenJS = js["grondlagen"].get<std::vector<json>>();
    for (unsigned int i = 0; i < grondlagenJS.size(); i++) {
        grondlagen.push_back(Grond(grondlagenJS[i]));
    }
    xPositie = js["x"].get<double>();
    yPositie = js["y"].get<double>();
    zettingsBerekeningJS = js;
    gen_msg();
}


Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType,float _xPos,float _yPos)
{
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = _yPos;
    gen_js();
    gen_msg();
}

void Zettingsberekening::gen_js(){
    json js;
    js["belasting"] = belastingsType.belastingsTypeJS;
    std::vector<json> grondlagen_js;
    for (unsigned int i = 0; i < grondlagen.size(); i++) {
        grondlagen_js.push_back(grondlagen[i].ground_js);
    }
    js["grondlagen"] = grondlagen_js;
    js["x"] = xPositie;
    js["y"] = yPositie;
    zettingsBerekeningJS = js;

}

Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType, float _xPos)
{
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = 0;
    gen_msg();

}

void Zettingsberekening::gen_msg()
{
    std::ostringstream str;
    str << "\nzettingsberekining in punt (" << std::setprecision(decimalPrecisionInShout)
        << xPositie << "," << std::setprecision(decimalPrecisionInShout) << yPositie << ")\n"
        << belastingsType.shout() << "\n";
    for (int i = 0; i < grondlagen.size(); i++) {
        str << grondlagen[i].shout() << "\n";
    }
    message = str.str();
    str.clear();
}

void Zettingsberekening::addGrondlaag(Grond g)
{
    grondlagen.push_back(g);
    message += g.shout()+"\n";
}

void Zettingsberekening::berekenZetting()
{
    dZettingPrim.clear();
    dSigma_eff.clear();
    dDelta_sigma.clear();
    graphDzetting.clear();
    double effectieveSpanningOpZ = 0;
    double diepte = 0;
    double totZetting = 0;
    double dDelt_sigma = 0;
    double finaleSpanningOpZ = 0, HOverC=0, lnGedeelte=0, zettingT=0;
    for (int i = 0; i < grondlagen.size(); i++) {
        double j = 0;
        double laagzetting = 0;
        while((j+(double)gridSize)<grondlagen[i].laagdikte){
            diepte += gridSize;
            effectieveSpanningOpZ+=grondlagen[i].drogeMassDichtheid*gridSize;
            dDelt_sigma = belastingsType.deltaSig(diepte, xPositie);
            finaleSpanningOpZ = effectieveSpanningOpZ + dDelt_sigma;
            HOverC = (gridSize) / (grondlagen[i].samendrukkingsCoeff);
            lnGedeelte = std::log( (finaleSpanningOpZ) / (double)effectieveSpanningOpZ);
            zettingT = (double)(HOverC*lnGedeelte);
            laagzetting += zettingT;
            //stockeer waarden voor latere preview
            dZettingPrim.push_back(zettingT);
            dDelta_sigma.push_back(dDelt_sigma);
            dSigma_eff.push_back(effectieveSpanningOpZ);
            j = j + gridSize;
        }
        grondlagen[i].primZetting = laagzetting;
    }
    double tot = 0;
    for (int k = 0; k < dZettingPrim.size(); k++) {
        tot += (double)dZettingPrim[k];
        graphDzetting.insert(graphDzetting.begin(), tot);
    }
    totalePrimaireZetting = tot;
}

void Zettingsberekening::berekenSecZetting(){
    //TODO: secundaire zetting, check params in grond. Als kan berekend worden-> bereken. anders? Zet dZettingSec in resp vector.
    for (int i = 0; i < grondlagen.size(); i++) {
        double t = 999999999; //neem t zeer groot, dit zal de secundaire zetting na een lange tijd berekenen, t_p moet kleiner zijn dan t
        grondlagen[i].secZetting= grondlagen[i].dikteNaPrim / (1 + grondlagen[i].e_p)*grondlagen[i].c_alpha*log(t / grondlagen[i].t_p);
    }

}

void Zettingsberekening::wijzigBelastingsType(BelastingsType b)
{
    belastingsType = b;
}

double Zettingsberekening::getTotaleZetting()
{
    return totalePrimaireZetting;
}

std::string Zettingsberekening::shout()
{
    std::ostringstream out;
    out<<message << std::endl;
    return out.str();
}

void Zettingsberekening::setGridSize(float _gridSize)
{
    gridSize = _gridSize;
}

float Zettingsberekening::getOpDiepte(float diepte, std::vector<double> &DiepteInfo)
{
    float Info_op_diepte = 0.0f;
    if (!DiepteInfo.empty()) {
        float n =(int)( (diepte / (float)gridSize));
        if (n >= 0 && n < DiepteInfo.size()) {
            Info_op_diepte = DiepteInfo[(int)n];
        }
    }
    else {
        std::cout << "diepte info was empty" << std::endl;
    }
    return Info_op_diepte;
}

float Zettingsberekening::getZettingOpDiepte(float diepte)
{
   
    return getOpDiepte(diepte,graphDzetting);
}

float Zettingsberekening::getEffectieveOpDiepte(float diepte)
{
    return getOpDiepte(diepte, dSigma_eff);
}

float Zettingsberekening::getDSigmaOpDiepte(float diepte)
{
    return getOpDiepte(diepte, dDelta_sigma);
}

float Zettingsberekening::getGridSize()
{
    return gridSize;
}

double Zettingsberekening::Consolidatiegraad(double Tv)
{
    double U = 0;
    bool uselessiteration = true;
    if (!uselessiteration) {
        if (Tv < 0.2827433389) {
            U = 2 * sqrt(PI*Tv) / PI;
        }
        else if (Tv >= 0.2827433389) {
            U = -0.01*exp(-2.467936863*Tv + 4.395395553) + 1;
        }
    }
    U = 1;
    double M = 0, expPart = 0, T2 = 0;
    for (double m = 0; m < sumPrecision; m++) {
        M = std::pow(0.5*PI*(2.0*m + 1), 2);
        expPart = std::exp(-M*Tv);
        std::cout << expPart << std::endl;
        T2 += (2 / M)*expPart;
    }
    U -= T2;
    return U;
}

double Zettingsberekening::Tijdsfactor(double U){
    double Tv = 0;
    if (U < 0.6) {
        Tv = PI / 4.0 * (U*U);
    }
    else {
        Tv = 1.781 - 0.933*std::log10(100 - 100 * U);
    }
    return Tv;
}

double Zettingsberekening::getZettingNaT(double t){
    //double T_v = c_v*t / (d*d);
    double bereikteZetting = 0;
    double TV = 0;
    for (unsigned int i = 0; i < grondlagen.size(); i++) {
        if (i > 0 && i < grondlagen.size() - 1) {
            //bepaal drainagelengte
            double drainagelength = getDrainageLength(grondlagen[i-1],
                grondlagen[i], grondlagen[i+1]);
            //bereken Tv = C_v*t/(D^2)
            TV = grondlagen[i].c_v*t / (drainagelength*drainagelength);
            //std::cout <<"T_v tussenin"<< TV << "\n";

            
        }
        else {
            //dus i =0 of i=grondlagen.size()-1
            //bereken direct T_v aangezien D=H/2 of dus D^2=H^2*0.25
            TV = grondlagen[i].c_v*t / (grondlagen[i].laagdikte*grondlagen[i].laagdikte*0.25);
            //std::cout << "T_v uiteinde" << TV << "\n";
        }
        //tel consolidatiegraad * totale primaire zetting vd laag = bereikte zetting vd laag
        bereikteZetting += Consolidatiegraad(TV)*grondlagen[i].primZetting;

    }
    return bereikteZetting;
}

double Zettingsberekening::getTimeToConsolidationDegree(double U){
    double ZettingBijU = U*totalePrimaireZetting;
    double TV = Tijdsfactor(U);
    return 0.0;
}

double Zettingsberekening::getDrainageLength(Grond & onder, Grond & huidig, Grond & boven){
    double onderKS =onder.k_s;
    double bovenKS = boven.k_s;
    double huidigKS = huidig.k_s;
    double factor = 1.0;
    if ( (bovenKS >= 100 * huidigKS) && (onderKS >= 100 * huidigKS) ){
        factor /= 2.0;
    }
    else if ((onderKS >= 100 * huidigKS) || (bovenKS >= 100 * huidigKS)
            && ! ((bovenKS >= 100 * huidigKS) && (onderKS >= 100 * huidigKS))  ) {
        factor /= 1.0;
    }

    return huidig.laagdikte*factor;
}


BelastingsType::BelastingsType()
{
    
}

BelastingsType::BelastingsType(json js){
    x1 = js["x1"].get<double>();
    x2 = js["x2"].get<double>();
    qs = js["qs"].get<double>();
    type = js["type"].get<int>();
    if (type == 0) {
        typeNaam = "Uniforme Stripbelasting";
    }

    belastingsBreedte = std::abs(x2 - x1);
    initialised = true;
    belastingsTypeJS = js;
}

BelastingsType::BelastingsType(float _x1,float _x2, float _qs, int _typeLast) :
    x1(_x1),x2(_x2)
    , qs(_qs), type(_typeLast)
{
    //breedte in m
    //qs is de blasting in kN/m^3
    //Type 0 : uniforme strip
    
    if (_typeLast == 0) {
        typeNaam = "Uniforme Stripbelasting";
    }
    
    belastingsBreedte = std::abs(x2 - x1);
    initialised = true;
    gen_js();

}

float BelastingsType::deltaSig(float z,float xPositie)
{
    float delta_sigma = 0;
    if (type == 0 && initialised) {
        //boussinesq uniforme stripbelasting 
        double gamma = std::atan((xPositie - x1) / z);
        double beta = std::atan((xPositie - x2) / z);
        double alpha = gamma - beta;
        delta_sigma = (qs / pi)*(alpha + sin(alpha)*cos(alpha+2*beta));
    }
    return delta_sigma;
}

void BelastingsType::gen_js(){
    json js;
    js["x1"] = x1;
    js["x2"] = x2;
    js["qs"] = qs;
    js["type"] = type;
    belastingsTypeJS = js;
}

std::string BelastingsType::shout()
{
    std::ostringstream out;
    
    out << "Belasting:\n " << std::setprecision(decimalPrecisionInShout)<<x1 << 
        "m tot " << std::setprecision(decimalPrecisionInShout)<< x2 << "m, \ngrootte "
        << std::setprecision(decimalPrecisionInShout)<<(qs) << "kN/m^2 \nType " << 
        typeNaam << std::endl;
    return out.str();
}
