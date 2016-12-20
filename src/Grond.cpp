#include <Grond.h>
Grond::Grond(float _samendrukkingscoeff, float _bovengrens, float _ondergrens,
             float _drogeMassaDichtheid, std::string _Naam, double c_v,
             double k_s = 0.0, double _nattemassadichtheid, double _OCR,
             double _ontlastingsconstante, double _c, double _c_a, double _phi,
             double _phi_a, double _c_alpha)
    : samendrukkingsCoeff(_samendrukkingscoeff),
      bovengrens(_bovengrens),
      Naam(_Naam),
      ondergrens(_ondergrens),
      drogeMassDichtheid(_drogeMassaDichtheid),
      c_v(c_v),
      k_s(k_s),
      natteMassadichtheid(_nattemassadichtheid),
      OCR(_OCR),
      ontlastingsconstante(_ontlastingsconstante),
      c(_c),
      c_a(_c_a),
      phi(_phi),
      phi_a(_phi_a),
      c_alpha(_c_alpha) {
    if (samendrukkingsCoeff == 0) {
        samendrukkingsCoeff = 0.0000000001;
    }
    if (natteMassadichtheid == 0) {
        natteMassadichtheid = drogeMassDichtheid * 1.2;
    }
    if (ontlastingsconstante == 1) {
        ontlastingsconstante = samendrukkingsCoeff;
    }
    laagdikte = std::abs(bovengrens - ondergrens);
    gen_js();
    gen_msg();
}

Grond::Grond(json grondJS) {
    samendrukkingsCoeff = grondJS["samendrukkingsCoeff"].get<double>();
    ondergrens = grondJS["ondergrens"].get<double>();
    bovengrens = grondJS["bovengrens"].get<double>();
    drogeMassDichtheid = grondJS["drogeMassaDichtheid"].get<double>();
    Naam = grondJS["Naam"].get<std::string>();
    natteMassadichtheid = grondJS["natteMassadichtheid"].get<double>();
    // natteMassadichtheid = drogeMassDichtheid;
    c_v = (grondJS["C_v"].get<double>() < 0.00000000000000001)
              ? 0
              : grondJS["C_v"].get<double>();
    k_s = (grondJS["k_s"].get<double>() < 0.0000000000000001)
              ? 0
              : grondJS["k_s"].get<double>();
    // protect against division by 0
    ontlastingsconstante =
        (grondJS["Ontlastingsconstante"].get<double>() < 0.00000000000000001)
            ? 1
            : grondJS["Ontlastingsconstante"].get<double>();
    OCR = (grondJS["OCR"].get<double>() < 1) ? 1 : grondJS["OCR"].get<double>();
    if (grondJS.find("c") != grondJS.end()) {
        c = grondJS["c"].get<double>();
    }
    if (grondJS.find("c_a") != grondJS.end()) {
        c_a = grondJS["c_a"].get<double>();
    }
    if (grondJS.find("phi") != grondJS.end()) {
        phi = grondJS["phi"].get<double>();
    }
    if (grondJS.find("phi_a") != grondJS.end()) {
        phi_a = grondJS["phi_a"].get<double>();
    }
    if (grondJS.find("c_alpha") != grondJS.end()) {
        c_alpha = grondJS["c_alpha"];
    }
    laagdikte = std::abs(bovengrens - ondergrens);
    // ground_js = grondJS;
    gen_js();
    gen_msg();
}

void Grond::gen_msg() {
    std::ostringstream str;

    str << "Grondlaag " << Naam << "\nsamendrukkingscoefficient C: "
        << std::setprecision(decimalPrecisionInShout) << samendrukkingsCoeff
        << "\nOntlastingsconstante: " << ontlastingsconstante << "\nOCR:" << OCR
        << " \n De laag gaat van " << std::setprecision(decimalPrecisionInShout)
        << (bovengrens) << "m tot "
        << std::setprecision(decimalPrecisionInShout) << (ondergrens)
        << "m\n droge massadichtheid: "
        << std::setprecision(decimalPrecisionInShout) << (drogeMassDichtheid)
        << "kN/m^3\nNatte massadichtheid: " << natteMassadichtheid
        << "kN/m^3\nmet C_v: " << c_v << "m^2/s\n"
        << "Doorlatendheid: " << k_s << "m/s\n"
        << "c: " << c << " phi: " << phi << " c': " << c_a
        << "  phi': " << phi_a << "\nc_alpha: " << c_alpha;
    str << std::endl;
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
    js["natteMassadichtheid"] = natteMassadichtheid;
    js["Naam"] = Naam;
    js["C_v"] = c_v;
    js["k_s"] = k_s;
    js["Ontlastingsconstante"] = ontlastingsconstante;
    js["OCR"] = OCR;
    js["c"] = c;
    js["phi"] = phi;
    js["phi_a"] = phi_a;
    js["c_a"] = c_a;
    js["c_alpha"] = c_alpha;
    ground_js = js;
}

Grond::~Grond() {}

std::string Grond::shout() { return Message; }

Zettingsberekening::~Zettingsberekening() {}

Zettingsberekening::Zettingsberekening() {}

Zettingsberekening::Zettingsberekening(json js) {
    belastingsType = BelastingsType(js["belasting"]);

    std::vector<json> grondlagenJS = js["grondlagen"].get<std::vector<json>>();
    for (unsigned int i = 0; i < grondlagenJS.size(); i++) {
        grondlagen.push_back(Grond(grondlagenJS[i]));
    }
    xPositie = js["x"].get<double>();
    yPositie = js["y"].get<double>();
    fea = js["fea"].get<double>();
    if (js.find("laagste FEA") != js.end()) {
        lowestPhea = js["laagste FEA"].get<double>();
    }
    gen_msg();
    gen_js();
    done = false;
}

Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType,
                                       float _xPos, float _yPos) {
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = _yPos;
    gen_js();
    gen_msg();
    done = false;
}

void Zettingsberekening::gen_js() {
    json js;
    js["belasting"] = belastingsType.belastingsTypeJS;
    std::vector<json> grondlagen_js;
    for (unsigned int i = 0; i < grondlagen.size(); i++) {
        grondlagen_js.push_back(grondlagen[i].ground_js);
    }
    js["grondlagen"] = grondlagen_js;
    js["x"] = xPositie;
    js["y"] = yPositie;
    js["fea"] = fea;
    js["laagste FEA"] = lowestPhea;
    zettingsBerekeningJS = js;
}

Zettingsberekening::Zettingsberekening(BelastingsType _belastingsType,
                                       float _xPos) {
    belastingsType = _belastingsType;
    xPositie = _xPos;
    yPositie = 0;
    gen_msg();
    done = false;
}

void Zettingsberekening::gen_msg() {
    std::ostringstream str;
    str << "\nzettingsberekining in punt ("
        << std::setprecision(decimalPrecisionInShout) << xPositie << ","
        << std::setprecision(decimalPrecisionInShout) << yPositie
        << ") Fea: " << fea << "m\n"
        << belastingsType.shout() << "\n";
    for (int i = 0; i < grondlagen.size(); i++) {
        str << grondlagen[i].shout() << "\n";
    }
    message = str.str();
    str.clear();
}

void Zettingsberekening::addGrondlaag(Grond g) {
    grondlagen.push_back(g);
    message += g.shout() + "\n";
    done = false;
}

void Zettingsberekening::berekenZetting() {
    if (!done) {
        dZettingPrim.clear();
        dSigma_eff.clear();
        dDelta_sigma.clear();
        graphDzetting.clear();
        int num_elements = 0;
        /*for (int i = 0; i < grondlagen.size(); ++i) {
            num_elements +=
                std::ceil(grondlagen[i].laagdikte / (double)gridSize);
            num_elements++;
        }*/
        num_elements = std::ceil(
            (grondlagen[0].bovengrens - grondlagen.back().ondergrens) /
            (double)gridSize);
        dZettingPrim.reserve(num_elements);
        dSigma_eff.reserve(num_elements);
        double num_elements2 =
            num_elements -
            (grondlagen.front().bovengrens - belastingsType.aanzetshoogte) /
                ((double)gridSize);
        dDelta_sigma.reserve(num_elements);
        double totaleSpanningOpZ = 0;
        double sigma_eff_op_aanzet = 0;
        bool corrected = false;

        maxZettingT = 0.00000001;
        double effectieveSpanningOpZ = 0;
        double diepte = 0;
        double totZetting = 0;
        double dDelt_sigma = 0;
        double critcalMassWeight = 0;
        double finaleSpanningOpZ = 0, HOverC = 0, lnGedeelte = 0, zettingT = 0;
        // placeholder for c and phi, find the lowest s_u
        double ESA_c_minimal = 0, ESA_phi_minimal = 0.00001,
               TSA_phi_minimal = 0.0001, s_u_TSA = 999999999999999,
               s_u_ESA = 9999999999999, TSA_c_minimal = 0;
        q_u_ESA = 999999;
        q_u_TSA = 999999;
        for (int i = 0; i < grondlagen.size(); i++) {
            double j = 0;
            double laagzetting = 0;
            // TODO: find the layer with the lowest s_u , this is the weakest
            // this is at least from the aanzetshoogte, downward. until a
            // certain depth has been reached. needs to check this
            // remove false when you want this to be activated
            // search till depth : aanzetshoogte - ??? is reached
            if (grondlagen[i].ondergrens < belastingsType.aanzetshoogte &&
                true) {
                // this might be wrong, the way s_u gets calculated p.e.
                // TSA su
                if (grondlagen[i].ondergrens > fea) {
                    critcalMassWeight = grondlagen[i].drogeMassDichtheid;
                } else {
                    critcalMassWeight =
                        grondlagen[i].natteMassadichtheid - waterGewicht;
                }
            }
            while ((j + (double)gridSize) < grondlagen[i].laagdikte) {
                diepte += gridSize;
                // als onder water
                if ((grondlagen[i].bovengrens - j) < fea) {
                    effectieveSpanningOpZ +=
                        (grondlagen[i].natteMassadichtheid - waterGewicht) *
                        gridSize;
                } else {
                    effectieveSpanningOpZ +=
                        (grondlagen[i].drogeMassDichtheid) * gridSize;
                }
                // find the total stresses with the lowest phea
                if ((grondlagen[i].bovengrens - j) < lowestPhea) {
                    totaleSpanningOpZ +=
                        (grondlagen[i].natteMassadichtheid) * gridSize;
                } else {
                    totaleSpanningOpZ +=
                        (grondlagen[i].drogeMassDichtheid) * gridSize;
                }
                // TODO check dit
                if ((grondlagen.front().bovengrens - diepte) <=
                    belastingsType.aanzetshoogte) {
                    if (!corrected) {
                        sigma_eff_op_aanzet = effectieveSpanningOpZ;
                        corrected = true;
                    }
                    dDelt_sigma =
                        belastingsType.deltaSig(
                            diepte - (((grondlagen.front().bovengrens -
                                        belastingsType.aanzetshoogte)) < 0
                                          ? 0
                                          : grondlagen.front().bovengrens -
                                                belastingsType.aanzetshoogte),
                            xPositie, yPositie) *
                        (belastingsType.qs - sigma_eff_op_aanzet);
                    dDelta_sigma.push_back(dDelt_sigma);
                }
                // hieronder. zal altijd false zijn als ocr =1 en ddelt_sigma >0
                if ((dDelt_sigma + effectieveSpanningOpZ) <
                    effectieveSpanningOpZ * grondlagen[i].OCR) {
                    HOverC = (gridSize) / (grondlagen[i].ontlastingsconstante);
                    lnGedeelte =
                        std::log((dDelt_sigma + effectieveSpanningOpZ) /
                                 (effectieveSpanningOpZ));
                    zettingT = (double)(HOverC * lnGedeelte);

                } else {
                    // deel belasting if not overgeconsolideerd
                    HOverC = (gridSize) / (grondlagen[i].samendrukkingsCoeff);
                    lnGedeelte =
                        std::log((dDelt_sigma + effectieveSpanningOpZ) /
                                 (effectieveSpanningOpZ * grondlagen[i].OCR));

                    zettingT = (double)(HOverC * lnGedeelte) +
                               (gridSize / grondlagen[i].ontlastingsconstante) *
                                   std::log(grondlagen[i].OCR);
                }
                finaleSpanningOpZ = effectieveSpanningOpZ + dDelt_sigma;
                // in ln : verhoogde op heersende
                laagzetting += zettingT;
                // stockeer waarden voor latere preview
                if (zettingT != 0) {
                    maxZettingT =
                        (zettingT > maxZettingT) ? zettingT : maxZettingT;
                    dZettingPrim.push_back(zettingT);
                }
                dSigma_eff.push_back(effectieveSpanningOpZ);
                dSigma_tot.push_back(totaleSpanningOpZ);
                j = j + gridSize;
            }
            grondlagen[i].primZetting = laagzetting;
            // std::cout<<laagzetting<<std::endl;
            if (grondlagen[i].ondergrens < belastingsType.aanzetshoogte &&
                (diepte - grondlagen[i].bovengrens +
                     belastingsType.aanzetshoogte <
                 2 * belastingsType.belastingsBreedte)) {
                double tempQUESA = calculateq_u(
                    grondlagen[i].c_a, grondlagen[i].phi_a, critcalMassWeight);
                if (tempQUESA < q_u_ESA && tempQUESA != 0) {
                    q_u_ESA = tempQUESA;
                }
            }
            if (grondlagen[i].ondergrens < belastingsType.aanzetshoogte) {
                double tempQUTSA = calculateq_u(
                    grondlagen[i].c, grondlagen[i].phi, critcalMassWeight);
                if (tempQUTSA < q_u_TSA && tempQUTSA != 0) {
                    q_u_TSA = tempQUTSA;
                }
            }
        }
        //        berekenSecZetting();

        // q_u_ESA =
        //     calculateq_u(ESA_c_minimal, ESA_phi_minimal, critcalMassWeight);
        double tot = 0;
        graphDzetting.resize(dZettingPrim.size(), 0.0);
        for (int k = 0; k < dZettingPrim.size(); k++) {
            tot += (double)dZettingPrim[k];
            graphDzetting[dZettingPrim.size() - 1 - k] = tot;
        }
        totalePrimaireZetting = tot;
        done = true;
        // writeToCSV();
    }
}

void Zettingsberekening::berekenSecZetting() {
    // TODO: secundaire zetting, check params in grond. Als kan berekend
    // worden-> bereken. anders? Zet dZettingSec in resp vector.
    for (int i = 0; i < grondlagen.size(); i++) {
        double t =
            9999999999999999;  // neem t zeer groot, dit zal de secundaire
        // zetting na een lange tijd berekenen, t_p moet
        // kleiner zijn dan t
        // find t_p
        grondlagen[i].getLaagDikteNaPrim();
        grondlagen[i].secZetting =
            grondlagen[i].dikteNaPrim / (1 + grondlagen[i].e_p) *
            grondlagen[i].c_alpha * log(t / grondlagen[i].t_p);
        tot_sec_zetting += grondlagen[i].secZetting;
    }
}

void Zettingsberekening::wijzigBelastingsType(BelastingsType b) {
    belastingsType = b;
    gen_msg();
    done = false;
}

void Zettingsberekening::writeToCSV() {
    if (!grondlagen.empty()) {
        std::string fileName = grondlagen.front().Naam;
        fileName += ".csv";
        std::ofstream targetFile;
        targetFile.open(fileName);
        targetFile << "sigmsEff,deltaSig,deltaZetting\n";
        for (unsigned int i = 0; i < dSigma_eff.size(); i++) {
            targetFile << dSigma_eff[i] << ",";
            int j = i - (dSigma_eff.size() - dDelta_sigma.size());
            if (j < 0) {
                targetFile << 0 << ",";
            } else {
                if (j < dDelta_sigma.size()) {
                    targetFile << dDelta_sigma[j] << ",";
                }
            }
            int k = i - (dSigma_eff.size() - dZettingPrim.size());
            if (k < 0) {
                targetFile << 0 << "\n";
            } else {
                if (k < dZettingPrim.size()) {
                    targetFile << dZettingPrim[k] << "\n";
                }
            }
        }
        targetFile.close();
    }
}

double Zettingsberekening::getTotaleZetting() { return totalePrimaireZetting; }

void Zettingsberekening::setPhea(double _fea) {
    fea = _fea;
    gen_msg();
    gen_js();
}

std::string Zettingsberekening::shout() {
    std::ostringstream out;
    out << message << std::endl;
    return out.str();
}

void Zettingsberekening::setGridSize(float _gridSize) {
    if (_gridSize != gridSize) {
        gridSize = _gridSize;
        done = false;
    }
}

float Zettingsberekening::getOpDiepte(float diepte,
                                      std::vector<double> &DiepteInfo) {
    float Info_op_diepte = 0.0f;
    if (!DiepteInfo.empty()) {
        float n = (int)((diepte / (float)gridSize));
        if (n >= 0 && n < DiepteInfo.size()) {
            Info_op_diepte = DiepteInfo[(int)n];
        }
    } else {
        std::cout << "diepte info was empty" << std::endl;
    }
    return Info_op_diepte;
}

float Zettingsberekening::getZettingOpDiepte(float diepte) {
    return getOpDiepte(diepte, graphDzetting);
}

float Zettingsberekening::getEffectieveOpDiepte(float diepte) {
    return getOpDiepte(diepte, dSigma_eff);
}

float Zettingsberekening::getDSigmaOpDiepte(float diepte) {
    float bovenGrensEersteLaag = 0;
    if (!grondlagen.empty()) {
        bovenGrensEersteLaag = grondlagen.front().bovengrens;
    }
    return getOpDiepte(
        (diepte - bovenGrensEersteLaag + belastingsType.aanzetshoogte),
        dDelta_sigma);
}

float Zettingsberekening::getTotSpanningOpDiepte(float diepte) {
    return getOpDiepte(diepte, dSigma_tot);
}

float Zettingsberekening::getGridSize() { return gridSize; }

double Zettingsberekening::Consolidatiegraad(double Tv) {
    double U = 0;
    bool iterationMethod = true;
    if (!iterationMethod) {
        if (Tv < 0.2827433389) {
            U = 2 * sqrt(PI * Tv) / PI;
        } else if (Tv >= 0.2827433389) {
            U = -0.01 * exp(-2.467936863 * Tv + 4.395395553) + 1;
        }
    }
    if (iterationMethod) {
        U = 1;
        double M = 0, expPart = 0, T2 = 0;
        double expTimesM = 1;
        for (double m = 0; m < sumPrecision && expTimesM != 0; m++) {
            M = std::pow(0.5 * PI * (2.0 * m + 1), 2);
            expPart = std::exp(-M * Tv);
            expTimesM = expPart * (2 / M);
            if (expTimesM == 0) {
                // std::cout << "aborting iteration after " << m << "
                // iterations"
                //           << std::endl;
            }
            T2 += expTimesM;
        }
        U -= T2;
    }
    return U;
}

void Zettingsberekening::setLowestPhea(double _phe) {
    lowestPhea = _phe;
    gen_msg();
    gen_js();
}

double Zettingsberekening::Tijdsfactor(double U) {
    double Tv = 0;
    if (U < 0.6) {
        Tv = PI / 4.0 * (U * U);
    } else {
        Tv = 1.781 - 0.933 * std::log10(100 - 100 * U);
    }
    return Tv;
}

double Zettingsberekening::getZettingNaT(double t) {
    // double T_v = c_v*t / (d*d);
    double bereikteZetting = 0;
    double TV = 0;
    for (unsigned int i = 0; i < grondlagen.size(); i++) {
        if (i > 0 && i < grondlagen.size() - 1) {
            // bepaal drainagelengte
            double drainagelength = getDrainageLength(
                grondlagen[i - 1], grondlagen[i], grondlagen[i + 1]);
            // bereken Tv = C_v*t/(D^2)
            TV = grondlagen[i].c_v * t / (drainagelength * drainagelength);
            // std::cout <<"T_v tussenin"<< TV << "\n";
        } else {
            // dus i =0 of i=grondlagen.size()-1
            // bereken direct T_v aangezien D=H/2 of dus D^2=H^2*0.25
            if (grondlagen.size() == 1) {
                // only 1 layer present
                // assume worst case: only seepage to the top
                TV = grondlagen[i].c_v * t / (pow(grondlagen[i].laagdikte, 2));
            } else {
                if (i == 0) {
                    // if the first layer
                    if (grondlagen[i].k_s <= 100 * grondlagen[i + 1].k_s) {
                        // seepage through layer below
                        TV = grondlagen[i].c_v * t /
                             pow(grondlagen[i].laagdikte * 0.5, 2);
                    } else {
                        TV = grondlagen[i].c_v * t /
                             pow(grondlagen[i].laagdikte, 2);
                    }
                } else {
                    // only the last one remains
                    if (grondlagen[i - 1].k_s >= 100 * grondlagen[i].k_s) {
                        TV = grondlagen[i].c_v * t /
                             pow(grondlagen[i].laagdikte, 2);
                    } else {
                        TV = grondlagen[i].c_v * t /
                             pow(grondlagen[i].laagdikte, 2);
                    }
                }
            }
            TV = grondlagen[i].c_v * t /
                 (grondlagen[i].laagdikte * grondlagen[i].laagdikte);
            // std::cout << "T_v uiteinde" << TV << "\n";
        }
        // tel consolidatiegraad * totale primaire zetting vd laag = bereikte
        // zetting vd laag
        bereikteZetting += Consolidatiegraad(TV) * grondlagen[i].primZetting;
    }
    return bereikteZetting;
}

double Zettingsberekening::getTimeToConsolidationDegree(double U) {
    double ZettingBijU = U * totalePrimaireZetting;
    double TV = Tijdsfactor(U);
    return 0.0;
}

double Zettingsberekening::getDrainageLength(Grond &onder, Grond &huidig,
                                             Grond &boven) {
    double onderKS = onder.k_s;
    double bovenKS = boven.k_s;
    double huidigKS = huidig.k_s;
    double factor = 1.0;
    if ((bovenKS >= 100 * huidigKS) && (onderKS >= 100 * huidigKS)) {
        factor /= 2.0;
    } else if ((onderKS >= 100 * huidigKS) ||
               (bovenKS >= 100 * huidigKS) &&
                   !((bovenKS >= 100 * huidigKS) &&
                     (onderKS >= 100 * huidigKS))) {
        factor /= 1.0;
    }

    return huidig.laagdikte * factor;
}

void Zettingsberekening::setPosition(double xCons, double yCons) {
    xPositie = (xCons);
    yPositie = (yCons);
    done = false;
    gen_msg();
}

BelastingsType::BelastingsType() {}

void BelastingsType::setAanzetshoogte(double z) {
    aanzetshoogte = z;
    gen_js();
}

BelastingsType::BelastingsType(json js) {
    x1 = js["x1"].get<double>();
    x2 = js["x2"].get<double>();
    qs = js["qs"].get<double>();
    type = js["type"].get<int>();
    if (type == 0) {
        typeNaam = "Uniforme plaat last";
    } else if (type == 1) {
        typeNaam = "Uniforme strip last";
    } else if (type == 2) {
        typeNaam = "uniform circular load";
        r = sqrt(x1 * x1 + x2 * x2);
    }
    if (js.find("aanzetshoogte") != js.end()) {
        aanzetshoogte = js["aanzetshoogte"].get<double>();
    }

    belastingsBreedte = std::abs(x2 - x1);
    initialised = true;
    belastingsTypeJS = js;
}

BelastingsType::BelastingsType(float _x1, float _x2, float _qs, int _typeLast,
                               double _aanzet)
    : x1(_x1), x2(_x2), qs(_qs), type(_typeLast) {
    // breedte in m
    // qs is de blasting in kN/m^3
    // Type 0 : uniforme strip

    if (_typeLast == 0) {
        typeNaam = "uniforme plaat last";
    } else if (_typeLast == 1) {
        typeNaam = "uniforme strip belasting";
    } else if (_typeLast == 2) {
        r = sqrt(x1 * x1 + x2 * x2);
        typeNaam = "uniform circular load";
    }
    aanzetshoogte = _aanzet;
    belastingsBreedte = std::abs(x2 - x1);
    initialised = true;
    gen_js();
}

float BelastingsType::deltaSig(float z, float xPositie, float yPositie) {
    float delta_sigma = 0;
    if (type == 1 && initialised) {
        // boussinesq uniforme stripbelasting
        double gamma = std::atan((xPositie - x1) / z);
        double beta = std::atan((xPositie - x2) / z);
        double alpha = gamma - beta;
        delta_sigma = (1 / pi) * (alpha + sin(alpha) * cos(alpha + 2 * beta));
    }

    if (type == 0 && initialised) {
        // subdivide in cases
        // case 1 x1 = L, x2=B resp in x and y dir., spawns rectangle with
        // diagonal (0,0)->(x1=L,x2=B)
        if (0 <= xPositie && xPositie <= x1 && 0 <= yPositie &&
            yPositie <= x2) {
            // geval 1.1
            double sig = 0.0;
            sig += sigma_plate_load(abs(xPositie), abs(x2 - yPositie), z);
            // geval 1.2
            sig += sigma_plate_load(abs(x1 - xPositie), abs(x2 - yPositie), z);
            // geval 1.3
            sig += sigma_plate_load(abs(xPositie), abs(yPositie), z);
            // geval 1.4
            sig += sigma_plate_load(abs(x1 - xPositie), abs(yPositie), z);
            delta_sigma = sig;
            return delta_sigma;
        }
        // case 2, x xor y out of rect load area
        else if (!(0 <= xPositie && xPositie <= x1) && 0 <= yPositie &&
                 yPositie <= x2) {
            double sig = 0.0;
            if (xPositie >= 0) {
                sig += sigma_plate_load(abs(xPositie), abs(x2 - yPositie), z);
                sig += sigma_plate_load(abs(xPositie), abs(yPositie), z);
                sig -=
                    sigma_plate_load(abs(xPositie - x1), abs(x2 - yPositie), z);
                sig -= sigma_plate_load(abs(xPositie - x1), abs(yPositie), z);
            } else {
                sig +=
                    sigma_plate_load(abs(xPositie - x1), abs(x2 - yPositie), z);
                sig += sigma_plate_load(abs(xPositie - x1), abs(yPositie), z);
                sig -= sigma_plate_load(abs(xPositie), abs(x2 - yPositie), z);
                sig -= sigma_plate_load(abs(xPositie), abs(yPositie), z);
            }
            delta_sigma = sig;
            return sig;
        } else if ((0 <= xPositie && xPositie <= x1) &&
                   !(0 <= yPositie && yPositie <= x2)) {
            // Copy paste of the previous case and change x and y
            double sig = 0.0;
            double xP = yPositie;
            double yP = xPositie;
            if (xP >= 0) {
                sig += sigma_plate_load(abs(xP), abs(x1 - yP), z);
                sig += sigma_plate_load(abs(xP), abs(yP), z);
                sig -= sigma_plate_load(abs(xP - x2), abs(x1 - yP), z);
                sig -= sigma_plate_load(abs(xP - x2), abs(yP), z);
            } else {
                sig += sigma_plate_load(abs(xP - x2), abs(x1 - yP), z);
                sig += sigma_plate_load(abs(xP - x2), abs(yP), z);
                sig -= sigma_plate_load(abs(xP), abs(x1 - yP), z);
                sig -= sigma_plate_load(abs(xP), abs(yP), z);
            }
            delta_sigma = sig;
            return sig;
        }
        // case 3, x and y both out of rectangular load area
        else if (!(0 <= xPositie && xPositie <= x1) &&
                 !(0 <= yPositie && yPositie <= x2)) {
            double sig = 0.0;
            double xPos = xPositie;
            double yPos = yPositie;

            if (xPos < 0 && yPos < 0) {
                xPos = x1 - xPos;
                yPos = x2 - yPos;
            } else if (xPos > 0 && yPos < 0) {
                yPos = x2 - yPos;
            } else if (xPos < 0 && yPos > 0) {
                xPos = x1 - xPos;
            }
            sig += sigma_plate_load(abs(xPos), abs(yPos), z);
            sig -= sigma_plate_load(abs(xPos), abs(yPos - x2), z);
            sig -= sigma_plate_load(abs(xPos - x1), abs(yPos), z);
            sig += sigma_plate_load(abs(xPos - x1), abs(yPos - x2), z);
            delta_sigma = sig;
            return sig;
        }
    }

    if (type == 2 && initialised) {
        delta_sigma = sigma_circular_load(z, r);
    }
    return delta_sigma;
}

double BelastingsType::sigma_circular_load(double z, double r) {
    double I0 = 1 - pow((1 / (1 + pow((r / z), 2))), 1.5);
    return I0;
}

void BelastingsType::gen_js() {
    json js;
    js["x1"] = x1;
    js["x2"] = x2;
    js["qs"] = qs;
    js["type"] = type;
    js["aanzetshoogte"] = aanzetshoogte;
    belastingsTypeJS = js;
}

std::string BelastingsType::shout() {
    std::ostringstream out;
    if (type == 0) {
        out << "Belasting:\n diagonaal (0,0)->("
            << std::setprecision(decimalPrecisionInShout) << x1 << ","
            << std::setprecision(decimalPrecisionInShout) << x2 << ") [m]";
    } else if (type == 1) {
        out << "van" << x1 << "tot" << x2;
    } else if (type == 2) {
        out << "met belastingsstraal" << sqrt(pow(x1, 2) + pow(x2, 2));
    }

    out << "\ngrootte " << std::setprecision(decimalPrecisionInShout) << (qs)
        << "kN/m^2 \nType " << typeNaam
        << "\nAanzetshoogte[m]: " << aanzetshoogte << std::endl;
    return out.str();
}

double BelastingsType::sigma_plate_load(double L, double B, double z) {
    double R1 = pow(L * L + z * z, 0.5);
    double R2 = pow(B * B + z * z, 0.5);
    double R3 = pow(L * L + B * B + z * z, 0.5);
    double s =
        1 / (2 * pi) * (atan(L * B / (z * R3)) +
                        L * B * z / R3 * (1 / pow(R1, 2) + 1 / pow(R2, 2)));
    return s;
}

double Zettingsberekening::calculateq_u(double c, double _phi,
                                        double massaGew) {
    // q_u = d_q*s_q*N_q*p_t+d_c*s_c*N_c*c+s_g*N_g*g_k*B/2
    // d_q = 1
    // s_q = 1+B/L * sin(phi')
    // s_c = (s_q*N_q-1)/(N_q-1)
    // s_g=1-0.3*B/L
    // s_q=s_c=s_g als L/B >5
    // N_q = exp(pi*tan(phi))*tan(Pi/4+phi/2)
    // N_c = (N_q-1)/tan(phi)
    // N1 = 2*(N_q-1)*tan(phi)
    double phi_inRads = _phi / 180 * pi;

    double phi_partialSafetyF = 1;
    double effectiveCohesion_safetyF = 1;
    if (phi_inRads == 0) {
        // then c represents s_u -> undrained shear strength;
        effectiveCohesion_safetyF = 1.4;
    }
    // hardcoded option to use safety factors
    if (true) {
        phi_partialSafetyF = 1.25;
        effectiveCohesion_safetyF = 1.25;
        if (phi_inRads == 0) {
            // then c represents s_u -> undrained shear strength;
            effectiveCohesion_safetyF = 1.4;
        }
    }

    double L = 1, B = 1;
    if (belastingsType.type == 0) {
        L = belastingsType.x1;
        B = belastingsType.x2;
    } else if (belastingsType.type = 1) {
        B = belastingsType.belastingsBreedte;
        L = 6 * B;
    } else if (belastingsType.type == 2) {
        B = belastingsType.r;
        L = belastingsType.r;
    }
    double s_q = 1;
    double s_g = 1;
    double s_c = 1;
    double N_q = exp(pi * tan(phi_inRads) / phi_partialSafetyF) *
                 tan(pi / 4.0 + phi_inRads / 2) *
                 tan(pi / 4.0 + phi_inRads / 2);
    double N_c = 0;
    if (phi_inRads != 0) {
        N_c = (N_q - 1) / tan(phi_inRads) * phi_partialSafetyF;
    } else if (phi_inRads == 0) {
        N_c = pi + 2;
    }
    double N_g = 2 * (N_q - 1) * tan(phi_inRads) / phi_partialSafetyF;
    double p_t = getEffectieveOpDiepte(grondlagen.front().bovengrens -
                                       belastingsType.aanzetshoogte);
    double d_c = 1;
    double d_q = 1;
    if (B != 0) {
        if (L / B < 5) {
            s_q = 1 + B / L * sin(phi_inRads);
            s_c = (s_q * N_q - 1) / (N_q - 1);
            s_g = 1 - 0.3 * B / L;
        }
    }
    double g_k = massaGew;

    double _qu = d_q * s_q * N_q * p_t +
                 d_c * s_c * N_c * c / effectiveCohesion_safetyF +
                 s_g * N_g * g_k * B / 2.0;
    // std::cout << "q_u = d_q*s_q*N_q*p_t+d_c*s_c*N_c*c+s_g*N_g*g_k*B/2\n"
    //           << d_q << "*" << s_q << "*" << N_q << "*" << p_t << "+" << d_c
    //           << "*" << s_c << "*" << N_c << "*" << c << "/" <<
    //           effectiveCohesion_safetyF  << "+" << s_g << "*"
    //           << N_g << "*" << g_k << "*" << B << "/" << 2.0 << "=" << _qu
    //           << "\n"
    //           << std::endl;
    return _qu;
}

double Zettingsberekening::getSU(double c, double phi, double sigma) {
    return c + sigma * abs(std::tan(phi / (180.0) * pi));
}
