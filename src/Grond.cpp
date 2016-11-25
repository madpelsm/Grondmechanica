#include <Grond.h>


Grond::Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,float _drogeMassaDichtheid, std::string _Naam):
    samendrukkingsCoeff(_samendrukkingscoeff),bovengrens(_bovengrens), Naam(_Naam),
    ondergrens(_ondergrens),drogeMassDichtheid(_drogeMassaDichtheid){
    std::replace(Naam.begin(), Naam.end(), ' ', '_');
    laagdikte = std::abs(bovengrens - ondergrens);


    std::ostringstream str;
    
    str<<"Grondlaag "<<Naam<<" met samendrukkingscoefficient C= " <<
        std::setprecision(decimalPrecisionInShout)<<samendrukkingsCoeff << 
        " \n De laag gaat van " << std::setprecision(decimalPrecisionInShout)<<(bovengrens)
        << "m tot " << std::setprecision(decimalPrecisionInShout)<<(ondergrens) << 
        "m \n De massadichtheid van de grond bedraagt: " << std::setprecision(decimalPrecisionInShout)
        <<(drogeMassDichtheid) << "kN/m^3";
    Message = str.str();
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


Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType,float _xPos,float _yPos)
{
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = _yPos;
    std::ostringstream str; 
    
    str<<"\nzettingsberekining in punt (" << std::setprecision(decimalPrecisionInShout)
        <<xPositie << "," << std::setprecision(decimalPrecisionInShout)<<yPositie<< ")\n" 
        << belastingsType.shout()<<"\n";
    message = str.str();

}

Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType, float _xPos)
{
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = 0;
    message = "\nzettingsberekining in punt (" + std::to_string(xPositie) + "," +
        std::to_string(yPositie) + ")\n";

}

void Zettingsberekening::addGrondlaag(Grond g)
{
    grondlagen.push_back(g);
    message += g.shout()+"\n";
}

void Zettingsberekening::berekenZetting()
{
    dZetting.clear();
    dSigma_eff.clear();
    dDelta_sigma.clear();
    double effectieveSpanningOpZ = 0;
    double diepte = 0;
    double totZetting = 0;
    double dDelt_sigma = 0;
    double finaleSpanningOpZ = 0, HOverC=0, lnGedeelte=0, zettingT=0;
    for (int i = 0; i < grondlagen.size(); i++) {
        double j = 0;
        while((j+(double)gridSize)<grondlagen[i].laagdikte){
            diepte += gridSize;
            effectieveSpanningOpZ+=grondlagen[i].drogeMassDichtheid*gridSize;
            dDelt_sigma = belastingsType.deltaSig(diepte, xPositie);
            finaleSpanningOpZ = effectieveSpanningOpZ + dDelt_sigma;
            HOverC = (gridSize) / (grondlagen[i].samendrukkingsCoeff);
            lnGedeelte = std::log( (finaleSpanningOpZ) / (double)effectieveSpanningOpZ);
            zettingT = (double)(HOverC*lnGedeelte);
            //stockeer waarden voor latere preview
            dZetting.push_back(zettingT);
            dDelta_sigma.push_back(dDelt_sigma);
            dSigma_eff.push_back(effectieveSpanningOpZ);
            j = j + gridSize;
        }
    }
    double tot = 0;
    for (int k = 0; k < dZetting.size(); k++) {
        tot += (double)dZetting[k];
        graphDzetting.insert(graphDzetting.begin(), tot);
    }
    totalePrimaireZetting = tot;
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


BelastingsType::BelastingsType()
{
    
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

std::string BelastingsType::shout()
{
    std::ostringstream out;
    
    out << "Belasting:\n " << std::setprecision(decimalPrecisionInShout)<<x1 << 
        "m tot " << std::setprecision(decimalPrecisionInShout)<< x2 << "m, \ngrootte "
        << std::setprecision(decimalPrecisionInShout)<<(qs) << "kN/m^2 \nType " << 
        typeNaam << std::endl;
    return out.str();
}
