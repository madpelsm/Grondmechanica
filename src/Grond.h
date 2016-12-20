#pragma once
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <json.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

using json = nlohmann::json;

class Grond {
    int decimalPrecisionInShout = 3;

   public:
    json ground_js;
    double samendrukkingsCoeff, bovengrens, ondergrens, laagdikte,
        drogeMassDichtheid, dikteNaPrim, secZetting, primZetting,
        natteMassadichtheid, OCR = 1, ontlastingsconstante = 0;
    double c = 999, c_a = 999, phi = 0, phi_a = 0;
    double e_p = 0, c_alpha = 0, t_p = 1, c_v, k_s;
    std::string Naam;
    std::string Message;
    Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,
          float _drogeMassaDichtheid, std::string _Naam, double c_v, double k_s,
          double _natteMassadichtheid = 0, double _OCR = 1,
          double _ontlastingsconstante = 1, double c = 999, double c_a = 999,
          double phi = 0, double phi_a = 0, double _c_alpha = 0);
    Grond(json grondJS);
    ~Grond();
    std::string shout();
    void gen_js();
    void gen_msg();
    double getLaagDikteNaPrim();
};

class BelastingsType {
    bool initialised = false;
    std::string typeNaam;
    int decimalPrecisionInShout = 3;

   public:
    double aanzetshoogte = 0;
    json belastingsTypeJS;
    float sigma_fin, sigmaZ_0;
    float belastingsBreedte, qs, pi = 4 * atan(1), x1, x2;
    double r = 1;
    int type = 0;  // type 0 : uniforme strip

    BelastingsType();
    BelastingsType(json js);
    BelastingsType(float x1, float x2, float _qs, int _typeLast,
                   double _aanzet = 0);
    float deltaSig(float z, float x, float yPos = 0);
    void gen_js();
    std::string shout();
    double sigma_plate_load(double L, double B, double z);
    double sigma_circular_load(double z, double r);
    void setAanzetshoogte(double z);
};

class Zettingsberekening {
    float gridSize = 0.1 /*gridgrootte in meter*/;
    double totalePrimaireZetting = 0;
    std::string message;
    int decimalPrecisionInShout = 3;

   public:
    bool done = false;
    json zettingsBerekeningJS;
    std::vector<double> dZettingPrim, dDelta_sigma, dSigma_eff, graphDzetting,
        dSigma_tot;
    BelastingsType
        belastingsType;  // bepaal adhv dit de formule voor de belasting
    double fea = 0, sumPrecision = 1000, PI = 4 * atan(1), waterGewicht = 9.81,
           lowestPhea = 0;
    double pi = 4 * atan(1);
    double q_u_ESA = 0, q_u_TSA = 0,
           maxZettingT = 0.000001;  // evenwichtsdraagvermogen
    double tot_sec_zetting = 0;
    float xPositie, yPositie;
    std::vector<Grond> grondlagen;
    ~Zettingsberekening();
    Zettingsberekening();
    Zettingsberekening(json js);
    Zettingsberekening(BelastingsType belastingsType, float _xPositie);
    Zettingsberekening(BelastingsType belastingsType, float _xPositie,
                       float _yPositie);
    void gen_js();
    void gen_msg();
    void addGrondlaag(Grond g);
    void berekenZetting();
    void setLowestPhea(double Fea_2);
    void berekenSecZetting();
    void wijzigBelastingsType(BelastingsType b);
    void writeToCSV();
    double getTotaleZetting();
    void setPhea(double fea);
    std::string shout();
    void setGridSize(float _gridSize);
    float getOpDiepte(float diepte,
                      std::vector<double> &DiepteInfo);  // bevat een vector de
                                                         // info van een param
                                                         // op verschillende
                                                         // dieptes
    float getZettingOpDiepte(float diepte);
    float getEffectieveOpDiepte(float diepte);
    float getDSigmaOpDiepte(float diepte);
    float getTotSpanningOpDiepte(float diepte);
    float getGridSize();
    double Consolidatiegraad(double Tv);  // use Tv (=Cv*t/D^2)
    double Tijdsfactor(double U);
    double getZettingNaT(double t);                 // using real time [s]
    double getTimeToConsolidationDegree(double U);  // use fraction of U;
    double getDrainageLength(Grond &onder, Grond &huidig, Grond &boven);
    void setPosition(double xCons, double yCons);
    double calculateq_u(double c, double phi, double massaGew);
    double getSU(double c, double phi, double sigma);
};
