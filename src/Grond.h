#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
class Grond {
    int decimalPrecisionInShout = 2;
public:
    double samendrukkingsCoeff, bovengrens, ondergrens, laagdikte, drogeMassDichtheid;
    std::string Naam;
    std::string Message;
    Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,float _drogeMassaDichtheid, std::string _Naam);
    ~Grond();
    std::string shout();
};


class BelastingsType {
    bool initialised = false;
    std::string typeNaam;
    int decimalPrecisionInShout = 2;
public:
    float sigma_fin, sigmaZ_0;
    float belastingsBreedte, qs, pi = 3.14159265359, x1, x2;
    int type = 0;//type 0 : uniforme strip

    BelastingsType();
    BelastingsType(float x1,float x2,float _qs,int _typeLast);
    float deltaSig(float z,float x);
    std::string shout();
};


class Zettingsberekening {
    float gridSize = 0.1/*gridgrootte in meter*/;
    double totalePrimaireZetting=0;
    std::string message;
    int decimalPrecisionInShout = 2;
public:
    std::vector<double>dZetting;
    BelastingsType belastingsType;//bepaal adhv dit de formule voor de belasting
    double fea;
    float  xPositie, yPositie;
    std::vector<Grond> grondlagen;
    ~Zettingsberekening();
    Zettingsberekening();
    Zettingsberekening(BelastingsType belastingsType, float _xPositie);
    Zettingsberekening(BelastingsType belastingsType,float _xPositie,float _yPositie);
    void addGrondlaag(Grond g);
    void berekenZetting();
    void wijzigBelastingsType(BelastingsType b);
    double getTotaleZetting();
    std::string shout();
    void setGridSize(float _gridSize);
    float getGridSize();
};