
#include <nanogui/nanogui.h>
#include <iostream>
#include <Grond.h>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <regex>
#include <thread>
#include <stdio.h>
#include <ctime>

using namespace nanogui;

enum test_enum {
    UniformStripLoad = 0,
    PuntLast = 1,
    LijnLast = 2
};

enum InfoOpDiepteTypes {
    Zetting = 0,
    EffectieveSpanning=1,
    SpanningsVerschil=2

};


//variable declarations

float gridSize = 0.001, graphScale = 1, infoDiepte = 0;
int windowMargin = 2;
bool multithreaded = true;
double bovengrens, ondergrens, samendrukkingsconstante, drogeMassadichtheid,
xPos, yPos = 0, sonderingsnummer, feaHoogte,
beginPosLast, eindPosLast, lastGrootte,c_v=0.0,k_s=0.0,tijd=999999;
int belastingID = 0;//de locatie van de belasting in de vector met belastingstypes
int WIDTH = 1000, HEIGHT = 700, pointPrecisionInDialog = 5;
std::string grondnaam = "Grondsoort";
InfoOpDiepteTypes diepteInfoType = Zetting;
test_enum loadType = UniformStripLoad;
Color grond_kleur_graph(0.543f,0.269f,0.0742f, 1.f);
std::vector<Zettingsberekening> sonderingsPunt;
std::vector<BelastingsType> belastingstypes;

//function identifiers
void addSolidLayerToSondering(double bovengrens, double ondergrens, 
        double samendrukkingsconstante, double drogemassadichtheid, 
        std::string grondnaam, int sondering, Screen * screen,double c_v,
        double k_s);
void genereerLast(double beginPosLast, double eindPosLast, double 
        lastGrootte, test_enum enumval, Screen * screen);
void genereerSonderingsPunt(int sonderingsnummer, double xPos, double yPos, 
        double feaHoogte, test_enum enumval, Screen * screen);
std::string writeToFile(std::string _content);
std::string makeSaveFile(std::vector<Zettingsberekening> &e);
void readFromFile();
std::vector<std::string> split(const std::string &s, char delim);

int main(int argc, char * argv[]) {
    nanogui::init();
    Screen *screen = new Screen(Vector2i(WIDTH, HEIGHT), "Zettingsberekening",true,false,8,8,24,8,2);
    

    bool enabled = true;

    FormHelper * belastingGUI = new FormHelper(screen);
    ref<Window> window3 = belastingGUI->addWindow(Eigen::Vector2i(windowMargin, 
            windowMargin), "Belastingsinformatie");

    FormHelper *grondGUI = new FormHelper(screen);
    ref<Window> window2 = grondGUI->addWindow(Eigen::Vector2i(3 * WIDTH / 5, HEIGHT / 2.6), "Grond");


    FormHelper *sonderingGUI = new FormHelper(screen);
    ref<Window> window1 = sonderingGUI->addWindow(Eigen::Vector2i(3.2*WIDTH / 5, windowMargin), "Sonderingsinformatie");

    FormHelper *zettingsGUI = new FormHelper(screen);
    ref<Window> window4 = zettingsGUI->addWindow(Eigen::Vector2i(windowMargin, HEIGHT / 5 * 2), "Zettingen");


    sonderingGUI->addGroup("Locatie")->setTooltip("Plaats waar we de zetting zullen berekenen.");
    sonderingGUI->addVariable("Zettingspunt", sonderingsnummer)->setTooltip(
            "Dit is het punt waar we de zetting in zullen berekenen, dit kan op de plaats van een sondering of een boring zijn.");
    sonderingGUI->addVariable("X positie [m]", xPos)->setTooltip("De plaats in de richting van de breedte van de last");
    sonderingGUI->addVariable("Y Positie [m]", yPos);
    sonderingGUI->addVariable("Freatisch oppervlakhoogte [m]", feaHoogte)->setTooltip(
            "Indien dit gekend is, dit is niet noodzakelijk");
    sonderingGUI->addButton("Maak zettingspunt aan", [&screen]() {
        genereerSonderingsPunt(sonderingsnummer, xPos, yPos, feaHoogte, loadType, screen);

    });

    grondGUI->addGroup("Grondlaag");
    grondGUI->addVariable("Grond naam", grondnaam)->setTooltip("Dit is niet noodzakelijk, maar juist invullen maakt controle eenvoudiger.");
    grondGUI->addVariable("Bovengrens [m]", bovengrens)->setTooltip(
            "Dit kan het aanvangspeil of ten opzichte T.A.W. zijn, zolang het onderling klopt, zal dit geen probleem opleveren.");
    grondGUI->addVariable("Ondergrens [m]", ondergrens)->setTooltip(
            "Tot waar de grondlaag reikt.");

    grondGUI->addGroup("Mechanische parameters");
    grondGUI->addVariable("Samendrukkingsconstante C", samendrukkingsconstante)->setTooltip(
            "Geprefereerd de C uit de laboproeven.");
    grondGUI->addVariable("Droge massadichtheid kN/m^3", drogeMassadichtheid);
    grondGUI->addVariable("C_v", c_v);
    grondGUI->addVariable("Doorlatendheidsfactor", k_s);

    grondGUI->addGroup("Bevestig");
    grondGUI->addButton("Voeg grondlaag toe", [&screen]() {
        addSolidLayerToSondering(bovengrens, ondergrens,
            samendrukkingsconstante, drogeMassadichtheid,
            grondnaam, sonderingsnummer, screen,c_v,k_s);
    });

    belastingGUI->addGroup("Belasting");
    belastingGUI->addVariable("Beginpositie van de last", beginPosLast);
    belastingGUI->addVariable("Eindpositie van de last", eindPosLast);
    belastingGUI->addVariable("Grote van de last [kN/m^3]", lastGrootte);
    belastingGUI->addVariable("Kies belastingsgeval", loadType, enabled)->setItems(
        { "Uniforme strip belasting","PuntLast","Lijnlast" });
    belastingGUI->addButton("Maak last aan", [&screen]() {
        genereerLast(beginPosLast, eindPosLast, lastGrootte, loadType, screen);
    })->setTooltip("Klik om de last aan te maken");
    int wW = 400, wH = 200;

    Window *w = new Window(screen, "Graph");
    w->setPosition(Eigen::Vector2i(2 * WIDTH / 6, 2.5 * HEIGHT / 7));
    w->setFixedWidth(wW);
    w->setVisible(true);


    TabWidget * graphTabs = new TabWidget(w);
    Widget * zettingGraphWidg = graphTabs->createTab("Zetting");
    Graph *graph = new Graph(zettingGraphWidg, "Zettingsverloop");
    graph->setFooter("diepte");
    VectorXf &func = graph->values();
    graph->setFixedHeight(wH);
    graph->setTooltip("Zettingsverloop");
    graph->setForegroundColor(grond_kleur_graph);



    Widget * sigma_eff_widg = graphTabs->createTab("sigma_eff");
    Graph * sigma_eff_graph = sigma_eff_widg->add<Graph>("sigma_eff");
    sigma_eff_graph->setFixedHeight(wH);
    sigma_eff_graph->setTooltip("Zettingsverloop");
    sigma_eff_graph->setForegroundColor(grond_kleur_graph);
    sigma_eff_graph->setFooter("diepte");
    VectorXf &sigma_eff_val = sigma_eff_graph->values();

    Widget * dSigm_widg = graphTabs->createTab("deltaSigma");
    Graph * dSigm_graph = dSigm_widg->add<Graph>("deltaSigma");
    dSigm_graph->setFixedHeight(wH);
    dSigm_graph->setTooltip("Zettingsverloop");
    dSigm_graph->setForegroundColor(grond_kleur_graph);
    dSigm_graph->setFooter("diepte");
    VectorXf &dSigm_val = dSigm_graph->values();

    Widget * dZettingWidg = graphTabs->createTab("Afgeleide zetting");
    Graph * dZetting_graph = dZettingWidg->add<Graph>("Afgeleide zetting");
    dZetting_graph->setFixedHeight(wH);
    dZetting_graph->setTooltip("Zettingsverloop");
    dZetting_graph->setForegroundColor(grond_kleur_graph);
    dZetting_graph->setFooter("diepte");
    VectorXf &dZetting_val = dZetting_graph->values();

    graphTabs->setLayout(new GridLayout(Orientation::Horizontal, 1, Alignment::Middle, 0, 0));

    zettingGraphWidg->setLayout(new GridLayout(Orientation::Vertical, 1, Alignment::Fill, 0, 0));
    dSigm_widg->setLayout(new GridLayout(Orientation::Vertical, 1, Alignment::Fill, 0, 0));
    sigma_eff_widg->setLayout(new GridLayout(Orientation::Vertical, 1, Alignment::Fill, 0, 0));
    dZettingWidg->setLayout(new GridLayout(Orientation::Vertical, 1, Alignment::Fill, 0, 0));
    
    PopupButton *popupBtn;
    popupBtn = new PopupButton(w, "", 0);
    popupBtn->setBackgroundColor(grond_kleur_graph);
    popupBtn->setFontSize(16);
    popupBtn->setFixedSize(Vector2i(70, 30));
    popupBtn->setPosition(Eigen::Vector2i((wW - 70) / 2.0, wH - 30));
    Popup * popup = popupBtn->popup();
    popup = popupBtn->popup();
    popup->setLayout(new GridLayout(Orientation::Horizontal, 2, Alignment::Fill, 0, 0));

    ColorWheel *colorwheel = new ColorWheel(popup);
    colorwheel->setColor(popupBtn->backgroundColor());
    Button *colorBtn = new Button(popup, "Kies");
    colorBtn->setFixedSize(Vector2i(100, 25));
    Color c = colorwheel->color();
    colorBtn->setBackgroundColor(c);
    colorwheel->setCallback([colorBtn](const Color &value) {

        colorBtn->setBackgroundColor(value);

    });
    colorBtn->setChangeCallback([colorBtn, popupBtn, &graph,&sigma_eff_graph,&dSigm_graph,&dZetting_graph](bool pushed) {
        if (pushed) {
            popupBtn->setBackgroundColor(colorBtn->backgroundColor());
            popupBtn->setPushed(false);
            graph->setForegroundColor(
                popupBtn->backgroundColor());
            sigma_eff_graph->setForegroundColor(
                popupBtn->backgroundColor());
            dSigm_graph->setForegroundColor(
                popupBtn->backgroundColor());
            dZetting_graph->setForegroundColor(
                popupBtn->backgroundColor());
        }
    });

    w->setLayout(new GridLayout(Orientation::Horizontal, 1, Alignment::Fill, 0, 0));
    graphTabs->setActiveTab(0);

    zettingsGUI->addVariable("Maasbreedte [m]", gridSize)->setTooltip(
            "Kies de nauwkeurigheid waarmee gerekend moet worden, kleinere waarden zorgen wel voor een langere rekentijd.");
    zettingsGUI->addButton("Calculate", [&screen, &func, &graph,&sigma_eff_val,&dSigm_val,&dZetting_val]() {
        std::ostringstream str;
        sonderingsPunt.size();
        str << "Onvolledige invoer :(";
        if (!sonderingsPunt.empty()) {
            str = std::ostringstream();
            //multi threaded calculation
            std::clock_t timer;
            timer = std::clock();
            double dur;
            if (multithreaded) {

                std::vector<std::thread> calculationThread;
                for (int i = 0; i < sonderingsPunt.size(); i++) {
                    sonderingsPunt[i].setGridSize(gridSize);
                    calculationThread.push_back(std::thread(&Zettingsberekening::berekenZetting, &sonderingsPunt[i]));
                }
                for (int i = 0; i < calculationThread.size(); i++) {
                    calculationThread[i].join();
                }
                calculationThread.clear();
            }
            for (int i = 0; i < sonderingsPunt.size(); i++) {
                if (!multithreaded) {
                    sonderingsPunt[i].setGridSize(gridSize);
                    sonderingsPunt[i].berekenZetting();
                }
                str << "zetting in punt (" << 
                    std::setprecision(pointPrecisionInDialog) << sonderingsPunt[i].xPositie << 
                    "," << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[i].yPositie 
                    << ") : " << std::to_string(sonderingsPunt[i].getTotaleZetting()) + "m.\nBerekend met maaswijdte: " 
                    << sonderingsPunt[i].getGridSize() << "m\n\n";
            }

            dur = (std::clock() - timer) / ((double)CLOCKS_PER_SEC);
            str << "Berekend in " << dur << "s\n" << std::endl;
            func.resize(sonderingsPunt[0].dZettingPrim.size());
            sigma_eff_val.resize(sonderingsPunt[0].dZettingPrim.size());
            dSigm_val.resize(sonderingsPunt[0].dZettingPrim.size());
            dZetting_val.resize(sonderingsPunt[0].dZettingPrim.size());
            for (int i = 0; i < sonderingsPunt[0].dZettingPrim.size(); i++) {
                dZetting_val[i] = graphScale*sonderingsPunt[0].dZettingPrim[i] / sonderingsPunt[0].dZettingPrim[0];
                dSigm_val[i] = graphScale*sonderingsPunt[0].dDelta_sigma[i] / sonderingsPunt[0].dDelta_sigma[0];
                sigma_eff_val[i]= graphScale*sonderingsPunt[0].dSigma_eff[i]/ sonderingsPunt[0].dSigma_eff.back();
                func[i] = graphScale*sonderingsPunt[0].graphDzetting[i]/sonderingsPunt[0].graphDzetting[0];
            }
            str.clear();

        }
        MessageDialog * m = new MessageDialog(screen, MessageDialog::Type::Information, "De zettingen", 
            str.str(), "OK", "Cancel");
    })->setTooltip("Zorg dat je zeker correcte input leverde.");
    zettingsGUI->addButton("Save", [&screen]() {
        MessageDialog *m = new MessageDialog(screen, MessageDialog::Type::Information, "Write status",
            writeToFile(makeSaveFile(sonderingsPunt)), "ok");
       
    });
    zettingsGUI->addButton("Load", []() {
        readFromFile();
    });

    zettingsGUI->addButton("Toon huidige configuratie", [&screen]() {
        std::string infoOfLoaded = "";
        /*for (int i = 0; i < sonderingsPunt.size(); i++) {
            infoOfLoaded += sonderingsPunt[i].shout();
        }*/
        if (sonderingsnummer < sonderingsPunt.size()) {
            infoOfLoaded += sonderingsPunt[sonderingsnummer].shout();
        }

        MessageDialog *m = new MessageDialog(screen, MessageDialog::Type::Information, "huidige configuratie", infoOfLoaded, "OK");
        })->setTooltip("Beschrijving van de grondlagen.");

    zettingsGUI->addButton("Export", []() {
        std::ostringstream str;
        for (int i = 0; i < sonderingsPunt.size(); i++) {
            str << sonderingsPunt[i].shout();
            str << "\nzetting in punt (" << std::setprecision(pointPrecisionInDialog) << 
                sonderingsPunt[i].xPositie << "," << std::setprecision(pointPrecisionInDialog) << 
                sonderingsPunt[i].yPositie << ") : " << std::to_string(sonderingsPunt[i].getTotaleZetting()) 
                + "mm.\nBerekend met maaswijdte: " << sonderingsPunt[i].getGridSize()<<"m\n********************************\n";

        }
        writeToFile(str.str());
    });

    zettingsGUI->addVariable("Diepte", infoDiepte);
    zettingsGUI->addVariable("Tijd [s]", tijd);
    zettingsGUI->addButton("Geef info op diepte", [&screen]() {

        if (sonderingsnummer < sonderingsPunt.size()) {
            std::ostringstream message;
            message << "Informatie op diepte " << std::setprecision(pointPrecisionInDialog) << infoDiepte << "m" << std::endl;
            std::ostringstream str;
            double zettingNaT = sonderingsPunt[sonderingsnummer].getZettingNaT(tijd);
            str << "Effectieve spanning: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getEffectieveOpDiepte(infoDiepte) << "kPa\n" <<
                "Spanningsverschil: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getDSigmaOpDiepte(infoDiepte) << "kPa\n" <<
                "Zetting: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getZettingOpDiepte(infoDiepte) << "m\n" << 
                "na "<<tijd<<"s is reeds: "<< zettingNaT <<"m zetting bereikt\n"
                <<"U="<<(zettingNaT/ sonderingsPunt[sonderingsnummer].getTotaleZetting())<<std::endl;
           
            MessageDialog* m2 = new MessageDialog(screen, MessageDialog::Type::Information, message.str(), str.str(), "OK", "Cancel", false);
        }
    });



    screen->setLayout(new GridLayout(Orientation::Vertical, 2, Alignment::Middle, windowMargin, 0));
    screen->setVisible(true);
    window3->setParent(screen);
    screen->performLayout();

    nanogui::mainloop();


    nanogui::shutdown();
    return 0;
}

void addSolidLayerToSondering(double bovengrens, double ondergrens, double samendrukkingsconstante,
    double drogemassadichtheid, std::string grondnaam, int sondering, Screen * screen,double c_v,double k_s) {

    if ((sondering) < sonderingsPunt.size()) {
        sonderingsPunt[sondering].addGrondlaag(Grond(samendrukkingsconstante, bovengrens, ondergrens, drogemassadichtheid, grondnaam,c_v,k_s));
        std::string t = sonderingsPunt[sondering].grondlagen[sonderingsPunt[sondering].grondlagen.size() - 1].shout();
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Information, "Grondlaag toegevoegd aan zettingsberekingspunt", t, "OK", "Cancel", false);
    }
    else {
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Warning, "Het sonderingspunt (" + std::to_string(sondering) + ") dat werd ingegeven is nog niet aangemaakt", "Gelieve eerst een sonderingspunt in te geven waar er een zettingsberekening dient uitgevoerd te worden", "OK", "Cancel", false);

    }
}

void genereerLast(double beginPosLast, double eindPosLast, double lastGrootte, test_enum enumval, 
    Screen * screen) {
    if (enumval < belastingstypes.size()) {
        // 1 belasting per keer. indien een nieuwe uniforme striplast -> overschrijf de vorige
        belastingstypes[enumval] = (BelastingsType(beginPosLast, eindPosLast, lastGrootte, enumval));
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Information, "Belasting aangemaakt", belastingstypes[enumval].shout(), "OK", "Cancel", false);
    }
    else if (enumval != 0) {
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Warning, "Belasting aangemaakt", "Het gekozen belastingstype wordt nog niet ondersteund", "OK", "Cancel", false);

    }
    else{
        belastingstypes.push_back(BelastingsType(beginPosLast, eindPosLast, lastGrootte, enumval));
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Information, "Belasting aangemaakt", belastingstypes[enumval].shout(), "OK", "Cancel", false);

    }
}

void genereerSonderingsPunt(int sonderingsnummer, double xPos, double yPos, double feaHoogte, 
    test_enum enumval, Screen * screen) {

    if (enumval < belastingstypes.size()) {
        if (sonderingsnummer == sonderingsPunt.size() ) {
            //dan nieuwe aanmaken
            sonderingsPunt.push_back(Zettingsberekening(belastingstypes[enumval], xPos, yPos));
            sonderingsPunt[sonderingsnummer].fea = feaHoogte;
            MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Information, "Zettingsberekeningspunt aangemaakt", sonderingsPunt[sonderingsnummer].shout(), "OK", "Cancel", false);

        }
        else if ((sonderingsnummer) < sonderingsPunt.size()) {
            sonderingsPunt[sonderingsnummer] = (Zettingsberekening(belastingstypes[enumval], xPos, yPos));
            sonderingsPunt[sonderingsnummer].fea = feaHoogte;
            MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Information, "Zettingsberekeningspunt aangemaakt", sonderingsPunt[sonderingsnummer].shout(), "OK", "Cancel", false);
        }
        else {
            MessageDialog *m = new MessageDialog(screen, MessageDialog::Type::Warning, "Oeps!", "De zettingsberekeningspunten dienen in een opeenvolgende volgorde ingevoerd te worden. Eg. eerst 0, dan 1, dan 2, enz..");
        }
    }
    else {
        MessageDialog* m = new MessageDialog(screen, MessageDialog::Type::Warning, "er is nog geen belastingstype aangemaakt", "Gelieve eerst een belastingstype aan te maken", "OK", "Cancel", false);
    }
}
std::string writeToFile(std::string fileContent) {

    std::string T = file_dialog({ { "SC","Simple Consolidation" },{ "txt","Textdocument" } }, true);
    if (!T.empty()) {
        std::ofstream myfile;
        myfile.open(T);
        myfile << fileContent;
        myfile.close();
        return "Write to file succesful";
    }
    return "write failed";
}

std::string makeSaveFile(std::vector<Zettingsberekening> &e) {

    json saveJS;
    std::vector<json> Zettingspuntenjs;
    for (int i = 0; i < e.size(); i++) {
        e[i].gen_js();
        Zettingspuntenjs.push_back(e[i].zettingsBerekeningJS);
    }
    saveJS["zettingsberekingspunten"] = Zettingspuntenjs;

    std::string t = saveJS.dump(4);
    return t;
}

void readFromFile() {
    std::string fileLocation = file_dialog({ { "SC","Simple Consolidation" } }, false);
    std::string Content;
    json inputJS;
    if (!fileLocation.empty()) {
        std::ifstream myfile;
        myfile.open(fileLocation);
        myfile >> inputJS;
        myfile.close();
    }
    std::vector<json> inputVectors = inputJS["zettingsberekingspunten"].get<std::vector<json>>();
    for (int i = 0; i < inputVectors.size(); i++) {
        sonderingsPunt.push_back(Zettingsberekening(inputVectors[i]));
        belastingstypes.push_back(sonderingsPunt[i].belastingsType);
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(std::move(item));
    }
    return elems;
}