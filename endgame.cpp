#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "endgame.h"
#define ALFA 0.0000001
#define BETA 0.0000005

using namespace std;

/*  comando da usare per compilare
/usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -static -s -o endgame grader.cpp endgame.cpp

PER MAC
/usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -o endgame endgame.cpp grader.cpp

PER WIN
/usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -static -s -o endgame endgame.cpp
timeout.exe 5 ./endgame

FORMULE

v = Vmax - W((Vmax-Vmin)/C)         v: velocita attuale                     C: capacita massima del guanto                                                    W:peso attuale delle pietre
E(p,t) = G(p) - R*T(p,t)            G(p): energia delle pietre raccolte     T(p,t): tempo impiegato a visitare seguendo cammino t e raccogliendo pietre p     R: fornito

T(per un determinato tratto) = lunghezza del tratto/Velocita' per percorrerlo l/v

*/

//------------ VARIABILI DEL PROBLEMA ------------

double E, G, T;       //E: energia finale                                 G: energia delle pietre raccolte        T:tempo impiegato
int N, S;             //N: numero città                                   S: città di partenza
int M, C;             //M: n° pietre                                      C: capacità del guanto
double R, Vmin, Vmax; //R: consumo energia guanto per unità di tempo      Vmin: velocità minima                   Vmax: velocità max Avengers
double Vpond;
ifstream in("input.txt");
ofstream out("output.txt");
set<int>::iterator it; //iteratore per SET di int
double bestE;          //miglior energia trovata fin'ora
int tmpC;              //capacita' attuale rimasta
vector<int> tmpPercorso, tmpRaccolta;
vector<int> pietreSulPercorso; //città i-esima che pietra ha preso
int counterPercorso = 0;
//vector< vector<double> > tableCoeff;

//------------ STRUTTURE ------------

struct Pietra
{
    double interesse;
    int counter;
    int peso;
    int energia;
    int index; //indice della pietra nel vector delle pietre
};

typedef struct Nodo
{
    vector<int> vic;        // (vicino, costo)
    bool visited;           //raggiunta o meno
    int nPietre;            //pietre in questa citta'
    vector<int> pietreNodo; //insieme delle pietre
    //double interesse;
    //double peso;
    //set<int> pietreNodo;        //insieme delle pietre
    //double media;
} Nodo;

vector<Nodo> grafo;
vector<Pietra> pietre;

//------------ METODI ------------

double max(double x, int y)
{
    if (x > y)
        return x;
    return y;
}

void print_solution()
{
    //3 numeri: Energia finale, Energia delle pietre, Tempo impiegato
    out << scientific << setprecision(10) << E << " ";
    out << scientific << setprecision(10) << G << " ";
    out << scientific << setprecision(10) << T << endl;
    //lista delle pietre raccolte in ogni citta, -1 se nessuna separate da spazi        print(bestPercorso)
    for (int i = 0; i < M; i++)
    {
        out << tmpRaccolta[i] << " ";
    }
    out << endl;
    //lista in ordine delle citta visitate, prima e ultima = S                          print(bestRaccolta)
    for (int i = 0; i < N + 1; i++)
    {
        out << tmpPercorso[i] << " ";
    }
    out << endl
        << "***" << endl;
}

void load_data()
{
    //cout<<"Inizio calcolo dati"<<endl;
    in >> N >> S;                      // 1) N: numero citta' S: posizione della base degli avengers
    in >> M >> C >> R >> Vmin >> Vmax; // 2) M: numero di pietre, C: capacita' del guanto, R:consumo di energia[double], Vmin e Vmax[double]

    //tableCoeff.resize(N, vector<double>(N));

    if (C != 0)
    {
        Vpond = ((Vmax - Vmin) / (double)C); //controllare che C non sia zero
    }
    else
    {
        Vpond = 0;
    }

    //cout<<"Resize vettori"<<endl;

    grafo.resize(N);
    pietre.resize(M);
    pietreSulPercorso.resize(N + 1);
    tmpRaccolta.resize(M);
    tmpC = C;

    //cout<<"Inizializza grafo"<<endl;

    for (int i = 0; i < N; i++) //inizializzazione a -1, nessuna pietra presa
    {
        grafo[i].vic.resize(N);
        pietreSulPercorso[i] = -1;
        grafo[i].vic[i] = -1;
        //grafo[i].interesse = 0;
    }
    pietreSulPercorso[N] = -1;

    //cout<<"Inizializza pietre"<<endl;

    for (int i = 0; i < M; i++) // M) m: massa della pietra i-sima, e: energia della pietra i-sima [int], [int]
    {
        in >> pietre[i].peso >> pietre[i].energia;
        pietre[i].index = i;
        pietre[i].interesse = pietre[i].energia / (double)pietre[i].peso;
        tmpRaccolta[i] = -1;
    }

    //cout<<"Inserisco pietre..."<<endl;
    int la;
    int indexCity;
    for (int i = 0; i < M; i++) // 2M) disponibilita' pietre, due righe ciascuna
    {                           //          -[int] lunghezza della lista
                                //          -lista di id delle citta che ospitano la pietra, [int int int]

        in >> la;
        pietre[i].counter = la; // quanto è frequente la pietra,

        for (int j = 0; j < la; j++)
        {
            in >> indexCity;
            grafo[indexCity].pietreNodo.push_back(i);
        }
    }

    //cout<<"Calcolo costi..."<<endl;
    int costx;
    for (int i = 1; i < N; i++) // N-1) matrice di adiacenza delle citta, solo parte inferiore con costi
    {

        double media;
        for (int j = 0; j < i; j++)
        {
            in >> costx;
            grafo[i].vic[j] = costx;
            grafo[j].vic[i] = costx;
        }
    }
}

void get_results()
{
    //v = Vmax - W((Vmax-Vmin)/C)
    double Vatt = Vmax;
    int PesoAtt = 0;
    double Eatt = 0.0;
    double Tatt = 0.0;

    E = 0.0; //E = G-T
    G = 0.0;
    T = 0.0;

    //E G T - calcolo

    int cammino = N;

    reverse(tmpPercorso.begin(), tmpPercorso.end());
    reverse(pietreSulPercorso.begin(), pietreSulPercorso.end());

    pietreSulPercorso[0] = pietreSulPercorso[N];
    pietreSulPercorso[N] = -1;

    for (int i = 0; i < cammino; i++)
    {
        if (pietreSulPercorso[i] != -1)
        {
            //cout << "Ho preso la pietra: " << pietreSulPercorso[i] << " nel nodo: " << tmpPercorso[i] << endl;
            PesoAtt += pietre[pietreSulPercorso[i]].peso;
            Eatt += (double)pietre[pietreSulPercorso[i]].energia;
        }
        if (C != 0) // calcolo della' velocita' in uscita dal nodo dopo aver preso la pietra
        {
            Vatt = (double)Vmax - (((double)PesoAtt) * Vpond);
        }
        int vicino = tmpPercorso[i + 1];

        //cout.precision(17);
        //cout<<"Tatt: "<<Tatt<< " \t\t\tdistanza:"<<grafo[tmpPercorso[i]].vic[vicino].second<<" Vatt: "<<Vatt<<" Nodo: "<<tmpPercorso[i]<<" j: "<<vicino<<" vicino:"<<grafo[tmpPercorso[i]].vic[vicino].first<<" i: "<<i<<endl;
        Tatt += ((double)(grafo[tmpPercorso[i]].vic[vicino])) / Vatt; // T(per un determinato tratto) = lunghezza del tratto/Velocita' per percorrerlo l/v
    }
    G = (double)Eatt;
    T = Tatt;
    E = G - R * T;

    //checkPietre();
    //E(p,t) = G(p) - R*T(p,t)
}
void checkPietre()
{
    for (int i = 0; i < M; i++)
    {
        if (tmpRaccolta[i] != -1)
        {
            //if (grafo[tmpRaccolta[i]].pietreNodo.count(i) != 1) {
            if (find(grafo[tmpRaccolta[i]].pietreNodo.begin(), grafo[tmpRaccolta[i]].pietreNodo.end(), i) == grafo[tmpRaccolta[i]].pietreNodo.end())
            {
                cout << "---------------------------ERRORE---------------------------" << endl;
                cout << "Hai preso la pietra: " << i << " nella città: " << tmpRaccolta[i] << " perché essa ha: " << endl;

                return;
            } /* else {
                cout<<"La pietra "<<i<<" è nel nodo "<<tmpRaccolta[i]<<" ; sue pietre ";
                for (it = grafo[tmpRaccolta[i]].pietreNodo.begin(); it != grafo[tmpRaccolta[i]].pietreNodo.end(); it++)
                {
                    cout<<*it<<" , ";
                }
                cout<<endl;
                
            } */
        }
    }
    //cout<<"finito"<<endl;
}

void alterSolution()
{
    reverse(tmpPercorso.begin(), tmpPercorso.end());
    reverse(pietreSulPercorso.begin(), pietreSulPercorso.end());

    pietreSulPercorso[0] = pietreSulPercorso[N];
    pietreSulPercorso[N] = -1;

    double worstInteresse = INT32_MAX;
    int worstIndex;
    int nodoAtt;
    //int bestPietra;
    //double bestIndex;
    int pietraAtt;

    for (int i = 0; i < N; i++) //trova pietra peggiore
    {
        if (pietreSulPercorso[i] != -1 && pietre[pietreSulPercorso[i]].interesse < worstInteresse)
        {
            worstInteresse = pietre[pietreSulPercorso[i]].interesse;
            worstIndex = i;
        }
    }

    tmpRaccolta[pietreSulPercorso[worstIndex]] = -1;
    tmpC += pietre[pietreSulPercorso[worstIndex]].peso;
    pietreSulPercorso[worstIndex] = -1;

    for (int i = worstIndex + 1; i < N; i++) //prova a sostituire
    {
        nodoAtt = tmpPercorso[i];
        //considera pietre migliori
        if (pietreSulPercorso[i] == -1)
        {
            if (grafo[nodoAtt].pietreNodo.size() > 0 && pietreSulPercorso[i] == -1)
            {
                for (int j = 0; j < grafo[nodoAtt].pietreNodo.size(); j++)
                {
                    pietraAtt = grafo[nodoAtt].pietreNodo[j];
                    if (tmpRaccolta[pietraAtt] == -1 && pietre[pietraAtt].interesse > worstInteresse && pietre[pietraAtt].peso <= tmpC)
                    {
                        tmpRaccolta[pietraAtt] = nodoAtt;
                        pietreSulPercorso[i] = pietraAtt;
                        tmpC -= pietre[pietraAtt].peso;
                        i = N;
                        break;
                    }
                }
            }
        }
    }
    //checkPietre();
    get_results();
}

double calcolaCoeff(int dist, int index) //intero che indica distanza, intero che indica indice nel vettore dei nodi
{
    int nPietreUtili = 0;
    double EnerTot = 0.0;
    int massaTot = 0;
    double coeff = 0.0;

    if (tmpC > 0)
    {

        for (int i = 0; i < M; i++)
        {
            //if (tmpRaccolta[i] == -1 && grafo[index].pietreNodo.count(i) == 1 && tmpC >= pietre[i].peso)
            if (tmpRaccolta[i] == -1 && find(grafo[index].pietreNodo.begin(), grafo[index].pietreNodo.end(), i) != grafo[index].pietreNodo.end() && tmpC >= pietre[i].peso)
            {
                nPietreUtili++;
                EnerTot += pietre[i].interesse;
                massaTot += pietre[i].peso;
            }
        }

        if (nPietreUtili == 0)
        {
            coeff = (double)-1 * dist;
        }
        else
        {
            EnerTot = EnerTot / (double)nPietreUtili;
            coeff = EnerTot / (double)dist;
        }
    }
    else
    {
        coeff = ((double)1) / ((double)dist);
    }

    return coeff;
}

double calcolaCoeffSoloDistanza(int dist, int index)
{
    double coeff = 0.0;
    coeff = ((double)1) / ((double)dist);
    return coeff;
}

void process_solution_greedy()
{

    //Calcolo in itinere
    /* double Vatt = Vmax;
    int PesoAtt = 0;
    double Eatt = 0.0;
    double Tatt = 0.0;

    E = 0.0; //E = G-T
    G = 0.0;
    T = 0.0; */

    set<int> nodiToVisit;
    for (int i = 0; i < N; i++) //inserisco tutti i nodi da visitare, tranne la sorgente
    {
        /* if (i == S) //skippa la sorgente
        {
            continue;
        } */
        nodiToVisit.insert(i);
    }

    int nodeHere = S;

    while (!nodiToVisit.empty()) //finche ci sono nodi da visitare
    {
        // cout << "Entra nel while" << endl;
        tmpPercorso.push_back(nodeHere);
        grafo[nodeHere].visited = true;
        vector<double> coeff; //vettore dei coefficienti dei vicini del nodo
        coeff.resize(grafo[nodeHere].vic.size());
        //  cout << "while-1" << endl;
        if (grafo[nodeHere].pietreNodo.size() > 0 && tmpC > 0)
        {

            int bestPietra = -1;
            int maxRapp = INT32_MIN;
            //      cout << "while-2" << endl;
            /* for (it = grafo[nodeHere].pietreNodo.begin(); it != grafo[nodeHere].pietreNodo.end(); ++it) //itera sulle pietre del nodo, controlla che non sia già stata presa globalmente e prende la migliore
            {
                //      cout << "Entrato nel for bestPietra" << endl;
                int pietra = *it;
                //      cout << "valore it: " << *it << endl;
                if (tmpRaccolta[pietra] == -1 && pietre[pietra].peso <= tmpC && pietre[pietra].interesse > maxRapp)
                {
                    bestPietra = pietra;
                    maxRapp = pietre[pietra].interesse;
                }
            } */
            vector<int>::iterator it;

            for (it = grafo[nodeHere].pietreNodo.begin(); it != grafo[nodeHere].pietreNodo.end(); ++it) //itera sulle pietre del nodo, controlla che non sia già stata presa globalmente e prende la migliore
            {
                //      cout << "Entrato nel for bestPietra" << endl;
                int pietra = *it;
                //      cout << "valore it: " << *it << endl;
                if (tmpRaccolta[pietra] == -1 && pietre[pietra].peso <= tmpC && pietre[pietra].interesse > maxRapp)
                {
                    bestPietra = pietra;
                    maxRapp = pietre[pietra].interesse;
                }
            }

            if (bestPietra != -1)
            {
                tmpRaccolta[bestPietra] = nodeHere;
                pietreSulPercorso[counterPercorso] = bestPietra;
                tmpC -= pietre[bestPietra].peso;

                //aggiunta peso ed energia
                /* PesoAtt += pietre[bestPietra].peso;
                Eatt += (double)pietre[bestPietra].energia;
                Vatt = (double)Vmax - (((double)PesoAtt) * Vpond); */
            }
        }
        int massimo = nodeHere;
        int massimoRiserva;
        //    cout << "Fuori bestPietra" << endl;
        if (nodiToVisit.size() != 1)
        {

            for (int i = 0; i < grafo[nodeHere].vic.size(); i++) //calcola coefficienti dei vicini
            {
                int vicino = i;
                int distanza = grafo[nodeHere].vic[i];

                if (vicino == S || grafo[vicino].visited || vicino == nodeHere)
                {
                    coeff[i] = INT32_MIN;
                    //cout<<"int min: "<<INT32_MIN<<endl;
                    continue;
                }
                else
                {
                    //coeff[i] = calcolaCoeffSoloDistanza(distanza, vicino);//tableCoeff[nodeHere][i];//calcolaCoeffSoloDistanza(distanza, vicino);
                    coeff[i] = calcolaCoeffSoloDistanza(distanza, vicino); //tableCoeff[nodeHere][i];//calcolaCoeffSoloDistanza(distanza, vicino);
                    /* if(distanza <= grafo[nodeHere].media) 
                    {
                        massimo = i; 
                        break;

                    } else {
                        massimoRiserva = i;
                    } */
                }
            }
            /* if (massimo == nodeHere)
            {
                massimo = massimoRiserva;
            } */

            for (int i = 0; i < coeff.size(); i++) //trovo l'indice dei vicini (non del grafo) che ha il coefficiente massimo
            {
                if (coeff[i] > coeff[massimo] && !grafo[i].visited)
                {
                    massimo = i;
                }
            }
            //cout<<"Massimo "<<massimo<<endl;
        }
        else
        {
            massimo = S;
        }
        nodiToVisit.erase(nodeHere);
        counterPercorso++;

        //Tatt += ((double)(grafo[nodeHere].vic[massimo])) / Vatt;

        nodeHere = massimo;
    }

    //cout << "While finito" << endl;
    tmpPercorso.push_back(S);
    get_results();
    //occhio all'ultimo link

    /* G = (double)Eatt;
    T = Tatt;
    E = G - R * T; */
}

void print_data()
{
    cout << " N: " << N << ", S: " << S << ", M: " << M << ", C: " << C << ", R: " << R << ", Vmax: " << Vmax << ", Vmin: " << Vmin << endl;
    cout << " --- informazioni pietre --- " << endl;
    for (int i = 0; i < M; i++)
    {
        cout << "Energia" << i << ": " << pietre[i].energia << ", "
             << "Massa" << i << ": " << pietre[i].peso << ", "
             << "Index" << i << ": " << pietre[i].index << endl;
    }
    cout << " --- informazioni grafo --- " << endl;

    for (int i = 0; i < N; i++)
    {
        cout << "Nodo " << i << " raggiunge : ";
        for (int j = 0; j < grafo[i].vic.size(); j++)
        {
            cout << j << " con costo " << grafo[i].vic[j] << " - ";
        }
        cout << endl;

        cout << "Nodo " << i << " ha n=" << grafo[i].pietreNodo.size() << " pietre : ";
        /* for (it = grafo[i].pietreNodo.begin(); it != grafo[i].pietreNodo.end(); it++)
        {
            cout << *it << " - ";
        } */
        cout << endl;
    }
}

int main()
{
    load_data();
    //print_data();
    process_solution_greedy();

    /* cout.precision(10);
    cout << "E: " << E << " G: " << G << " T: " << T << " R:" << R << endl;

    cout << "Percorso: " << endl;
    for (int x : tmpPercorso)
    {
        cout << x << " | ";
    }
    cout << endl
         << "Raccolta: " << endl;

    for (int x : tmpRaccolta)
    {
        cout << x << " | ";
    }
    cout << endl; */
    bestE = E;
    print_solution();
    while (true)
    {
        alterSolution();
        if (E > bestE)
        {
            bestE = E;
            /* bestG = G;
            bestT = T;
            bestPercorso = tmpPercorso;
            bestRaccolta = tmpRaccolta; */
            print_solution();
        }
    }

    return 0;
}