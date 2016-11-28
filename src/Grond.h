#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <json.hpp>
#include <algorithm>

using json = nlohmann::json;

class Grond {
    int decimalPrecisionInShout = 3;
public:
    json ground_js;
    double samendrukkingsCoeff, bovengrens, ondergrens, laagdikte, drogeMassDichtheid;
    std::string Naam;
    std::string Message;
    Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,float _drogeMassaDichtheid, std::string _Naam);
    Grond(json grondJS);
    ~Grond();
    std::string shout();
    void gen_js();
    void gen_msg();
};


class BelastingsType {
    bool initialised = false;
    std::string typeNaam;
    int decimalPrecisionInShout = 3;
public:
    json belastingsTypeJS;
    float sigma_fin, sigmaZ_0;
    float belastingsBreedte, qs, pi = 4*atan(1), x1, x2;
    int type = 0;//type 0 : uniforme strip

    BelastingsType();
    BelastingsType(json js);
    BelastingsType(float x1,float x2,float _qs,int _typeLast);
    float deltaSig(float z,float x);
    void gen_js();
    std::string shout();
};


class Zettingsberekening {
    float gridSize = 0.1/*gridgrootte in meter*/;
    double totalePrimaireZetting=0;
    std::string message;
    int decimalPrecisionInShout = 3;
public:
    json zettingsBerekeningJS;
    std::vector<double>dZetting,dDelta_sigma,dSigma_eff,graphDzetting;
    BelastingsType belastingsType;//bepaal adhv dit de formule voor de belasting
    double fea;
    float  xPositie, yPositie;
    std::vector<Grond> grondlagen;
    ~Zettingsberekening();
    Zettingsberekening();
    Zettingsberekening(json js);
    Zettingsberekening(BelastingsType belastingsType, float _xPositie);
    Zettingsberekening(BelastingsType belastingsType,float _xPositie,float _yPositie);
    void gen_js();
    void gen_msg();
    void addGrondlaag(Grond g);
    void berekenZetting();
    void wijzigBelastingsType(BelastingsType b);
    double getTotaleZetting();
    std::string shout();
    void setGridSize(float _gridSize);
    float getOpDiepte(float diepte, std::vector<double>& DiepteInfo);//bevat een vector de info van een param op verschillende dieptes
    float getZettingOpDiepte(float diepte);
    float getEffectieveOpDiepte(float diepte);
    float getDSigmaOpDiepte(float diepte);
    float getGridSize();
};