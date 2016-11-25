
#include <nanogui/nanogui.h>
#include <iostream>
#include <Grond.h>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <regex>

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
int windowMargin = 8;
double bovengrens, ondergrens, samendrukkingsconstante, drogeMassadichtheid,
xPos, yPos = 0, sonderingsnummer, feaHoogte,
beginPosLast, eindPosLast, lastGrootte;
int belastingID = 0;//de locatie van de belasting in de vector met belastingstypes
int WIDTH = 800, HEIGHT = 600, pointPrecisionInDialog = 3;
std::string grondnaam = "Grondsoort";
InfoOpDiepteTypes diepteInfoType = Zetting;
test_enum loadType = UniformStripLoad;
Color grond_kleur_graph(0.543f,0.269f,0.0742f, 1.f);
std::vector<Zettingsberekening> sonderingsPunt;
std::vector<BelastingsType> belastingstypes;

//function identifiers
void addSolidLayerToSondering(double bovengrens, double ondergrens, 
        double samendrukkingsconstante, double drogemassadichtheid, 
        std::string grondnaam, int sondering, Screen * screen);
void genereerLast(double beginPosLast, double eindPosLast, double 
        lastGrootte, test_enum enumval, Screen * screen);
void genereerSonderingsPunt(int sonderingsnummer, double xPos, double yPos, 
        double feaHoogte, test_enum enumval, Screen * screen);
std::string writeToFile(std::string _content);
std::string makeSaveFile(std::vector<Zettingsberekening> e);
void readFromFile();
std::vector<std::string> split(const std::string &s, char delim);

int main(int argc, char * argv[]) {
    nanogui::init();
    Screen *screen = new Screen(Vector2i(WIDTH, HEIGHT), "Zettingsberekening");


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
    sonderingGUI->addVariable("X positie [m]", xPos);
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

    grondGUI->addGroup("Bevestig");
    grondGUI->addButton("Voeg grondlaag toe", [&screen]() {
        addSolidLayerToSondering(bovengrens, ondergrens,
            samendrukkingsconstante, drogeMassadichtheid,
            grondnaam, sonderingsnummer, screen);
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
    int wW = 200, wH = 300;
    Window *w = new Window(screen, "Graph");
    w->setPosition(Eigen::Vector2i(2 * WIDTH / 6, 3 * HEIGHT / 7));
    w->setWidth(wW);
    w->setHeight(wH);
    w->setVisible(true);
    Graph *graph = new Graph(w, "Zettingsverloop");
    graph->setFooter("diepte");
    VectorXf &func = graph->values();
    graph->setPosition(Eigen::Vector2i(0, 30));
    graph->setFixedSize(Eigen::Vector2i(wW, wH - 60));
    graph->setTooltip("Zettingsverloop");
    graph->setForegroundColor(grond_kleur_graph);
    PopupButton *popupBtn;

    popupBtn = new PopupButton(w, "", 0);

    popupBtn->setBackgroundColor(grond_kleur_graph);

    popupBtn->setFontSize(16);
    popupBtn->setFixedSize(Vector2i(70, 30));
    popupBtn->setPosition(Eigen::Vector2i((wW - 70) / 2.0, wH - 30));
    Popup * popup = popupBtn->popup();
    popup = popupBtn->popup();

    popup->setLayout(new GroupLayout());



    ColorWheel *colorwheel = new ColorWheel(popup);

    colorwheel->setColor(popupBtn->backgroundColor());

    Button *colorBtn = new Button(popup, "Kies");

    colorBtn->setFixedSize(Vector2i(100, 25));

    Color c = colorwheel->color();

    colorBtn->setBackgroundColor(c);



    colorwheel->setCallback([colorBtn](const Color &value) {

        colorBtn->setBackgroundColor(value);

    });
    colorBtn->setChangeCallback([colorBtn, popupBtn, &graph](bool pushed) {

        if (pushed) {

            popupBtn->setBackgroundColor(colorBtn->backgroundColor());

            popupBtn->setPushed(false);

            graph->setForegroundColor(
                popupBtn->backgroundColor());

        }

    });


    zettingsGUI->addVariable("Maasbreedte [m]", gridSize)->setTooltip(
            "Kies de nauwkeurigheid waarmee gerekend moet worden, kleinere waarden zorgen wel voor een langere rekentijd.");
    zettingsGUI->addButton("Calculate", [&screen, &func, &graph]() {
        std::ostringstream str;
        sonderingsPunt.size();
        str << "Onvolledige invoer :(";
        if (!sonderingsPunt.empty()) {
            str = std::ostringstream();
            for (int i = 0; i < sonderingsPunt.size(); i++) {
                sonderingsPunt[i].setGridSize(gridSize);
                sonderingsPunt[i].berekenZetting();
                str << "zetting in punt (" << 
                    std::setprecision(pointPrecisionInDialog) << sonderingsPunt[i].xPositie << 
                    "," << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[i].yPositie 
                    << ") : " << std::to_string(sonderingsPunt[i].getTotaleZetting()) + "m.\nBerekend met maaswijdte: " 
                    << sonderingsPunt[i].getGridSize() << "m\n";
            }
            func.resize(sonderingsPunt[0].grondlagen.size());
            for (int i = 0; i < sonderingsPunt[0].grondlagen.size(); i++) {
                if (i < 0) std::cout << "negative iterator" << std::endl;
                func[i] = graphScale*sonderingsPunt[0].dZetting[i]/sonderingsPunt[0].dZetting[0];
            }

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
        for (int i = 0; i < sonderingsPunt.size(); i++) {
            infoOfLoaded += sonderingsPunt[i].shout();
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

    //zettingsGUI->addVariable("Kies het infotype", infoOpDiepte, true)->setItems({"De zetting","Effectieve spanning","Spanningsverschil"});
    zettingsGUI->addVariable("Diepte", infoDiepte);
    zettingsGUI->addButton("Geef info op diepte", [&screen]() {

        if (sonderingsnummer < sonderingsPunt.size()) {
            std::ostringstream message;
            message << "Informatie op diepte " << std::setprecision(pointPrecisionInDialog) << infoDiepte << "m" << std::endl;
            std::ostringstream str;
            str << "Effectieve spanning: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getEffectieveOpDiepte(infoDiepte) << "kPa\n" <<
                "Spanningsverschil: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getDSigmaOpDiepte(infoDiepte) << "kPa\n" <<
                "Zetting: " << std::setprecision(pointPrecisionInDialog) << sonderingsPunt[sonderingsnummer].getZettingOpDiepte(infoDiepte) << "m" << std::endl;
           
            MessageDialog* m2 = new MessageDialog(screen, MessageDialog::Type::Information, message.str(), str.str(), "OK", "Cancel", false);
        }
    });
/*
    // tabs
    Window * window = new Window(screen, "titel");
    window->setLayout(new GroupLayout());
    window->setPosition(Vector2i(425, 15));
    window->setSize(Eigen::Vector2i(300, 500));


    TabWidget* tabWidget = window->add<TabWidget>();
    Widget* layer = tabWidget->createTab("Color WheelTAV");

    layer->setLayout(new GroupLayout());
    layer->add<Label>("Color wheel widget", "sans-bold");
    layer->add<ColorWheel>();

    Widget *layer2 = tabWidget->createTab("function graph");

    layer2->setLayout(new GroupLayout());
    layer2->add<Label>("Function graph widget", "sans-bold");
    graph->setParent(layer2);
     */
    screen->setVisible(true);
    window3->setParent(screen);
    screen->performLayout();

    nanogui::mainloop();


    nanogui::shutdown();
    return 0;
}

void addSolidLayerToSondering(double bovengrens, double ondergrens, double samendrukkingsconstante,
    double drogemassadichtheid, std::string grondnaam, int sondering, Screen * screen) {

    if ((sondering) < sonderingsPunt.size()) {
        sonderingsPunt[sondering].addGrondlaag(Grond(samendrukkingsconstante, bovengrens, ondergrens, drogemassadichtheid, grondnaam));
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

std::string makeSaveFile(std::vector<Zettingsberekening> e) {
    std::string t;
    for (int i = 0; i < e.size(); i++) {
        t += "x" + std::to_string(e[i].xPositie) + "y" + std::to_string(e[i].yPositie) + "w" + std::to_string(e[i].fea);
        t += "b" + std::to_string(e[i].belastingsType.type) + "bxone" + 
            std::to_string(e[i].belastingsType.x1) + "bxtw" + std::to_string(e[i].belastingsType.x2) + 
            "bq" + std::to_string(e[i].belastingsType.qs);
        for (int j = 0; j < e[i].grondlagen.size(); j++) {
            t += "gyup" + std::to_string(e[i].grondlagen[j].bovengrens) + "gyun" + 
                std::to_string(e[i].grondlagen[j].ondergrens) + "C" + std::to_string(e[i].grondlagen[j].samendrukkingsCoeff) + 
                "gm" + std::to_string(e[i].grondlagen[j].drogeMassDichtheid) + "gn" + e[i].grondlagen[j].Naam + "$";
        }
        t += ",";
    }
    return t;
}

void readFromFile() {
    std::string fileLocation = file_dialog({ { "SC","Simple Consolidation" } }, false);
    std::string Content;
    if (!fileLocation.empty()) {
        std::ifstream myfile;
        myfile.open(fileLocation);
        myfile >> Content;
        myfile.close();
    }
    //parse to usefull data
    std::vector<std::string> sonderingsPunten = split(Content, ',');
    std::vector<Zettingsberekening> importedZettingsPunten;
    std::vector<BelastingsType> importedBelasting;
    for (int i = 0; i < sonderingsPunten.size(); i++) {

        std::size_t foundX = sonderingsPunten[i].find("x");
        std::size_t foundY = sonderingsPunten[i].find("y", 0);
        std::size_t foundW = sonderingsPunten[i].find("w", 0);
        std::size_t foundB = sonderingsPunten[i].find("b", 0);
        std::size_t foundBXONE = sonderingsPunten[i].find("bxone", 0);
        std::size_t foundBXTW = sonderingsPunten[i].find("bxtw", 0);
        std::size_t foundBQ = sonderingsPunten[i].find("bq", 0);
        std::size_t foundGYUP = sonderingsPunten[i].find("gyup", 0);


        double x = std::stod(sonderingsPunten[i].substr(foundX + 1, (foundY - foundX) - 1));

        double y = std::stod(sonderingsPunten[i].substr(foundY + 1, (foundW - foundY) - 1));
        double loaded_fea = std::stod(sonderingsPunten[i].substr(foundW + 1, foundB - foundW - 1));
        double belastings_type = std::stod(sonderingsPunten[i].substr(foundB + 1, foundBXONE - foundB - 1));
        double belasting_x1 = std::stod(sonderingsPunten[i].substr(foundBXONE + 5, foundBXTW - foundBXONE - 5));
        double belasting_x2 = std::stod(sonderingsPunten[i].substr(foundBXTW + 4, foundBQ - foundBXTW - 4));
        double belasting_qs = std::stod(sonderingsPunten[i].substr(foundBQ + 2, foundGYUP - foundBQ - 2));
        importedBelasting.push_back(BelastingsType(belasting_x1, belasting_x2, belasting_qs, belastings_type));
        importedZettingsPunten.push_back(Zettingsberekening(BelastingsType(belasting_x1, belasting_x2, 
            belasting_qs, belastings_type), x, y));


        int lengte = sonderingsPunten[i].size() - foundGYUP;
        std::vector<std::string> grond_lagen = split(sonderingsPunten[i].substr(foundGYUP, lengte - 1), '$');
        for (int j = 0; j < grond_lagen.size(); j++) {
            std::size_t foundGyup = grond_lagen[j].find("gyup", 0);
            std::size_t foundGYUN = grond_lagen[j].find("gyun", 0);
            std::size_t foundC = grond_lagen[j].find("C", 0);
            std::size_t foundGM = grond_lagen[j].find("gm", 0);
            std::size_t foundGN = grond_lagen[j].find("gn", 0);

            double ground_yUp = std::stod(grond_lagen[j].substr(foundGyup + 4, foundGYUN - foundGyup - 4));
            double ground_yUnder = std::stod(grond_lagen[j].substr(foundGYUN + 4, foundC - foundGYUN - 4));
            double C = std::stod(grond_lagen[j].substr(foundC + 1, foundGM - foundC - 1));
            double dryMass = std::stod(grond_lagen[j].substr(foundGM + 2, foundGN - foundGM - 2));
            std::string grondNaam = grond_lagen[j].substr(foundGN + 2, grond_lagen[j].size() - foundGN);
            importedZettingsPunten[i].addGrondlaag(Grond(C, ground_yUp, ground_yUnder, dryMass, grondNaam));
        }
    }
    if (importedZettingsPunten.size() > 0) {
        sonderingsPunt = importedZettingsPunten;
    }
    if (importedBelasting.size() > 0) {
        belastingstypes = importedBelasting;
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