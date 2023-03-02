#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
//#include "endgame.h"

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
ifstream in("input11.txt");
ofstream out("output.txt");
set<int>::iterator it; //iteratore per SET di int

double bestE, bestG, bestT;
int tmpC; //capacita' attuale rimasta
vector<int> bestPercorso, bestRaccolta, tmpPercorso, tmpRaccolta;

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
    vector<pair<int, int> > vic; // (vicino, costo)
    bool visited;               //raggiunta o meno
    int nPietre;                //pietre in questa citta'
    set<int> pietreNodo;        //insieme delle pietre
} Nodo;

vector<Nodo> grafo;
vector<Pietra> pietre;
vector<int> pietreInteressanti;

//------------ METODI ------------

int max(int x, int y)
{
    return (x > y) ? x : y;
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


bool sortPietre(int p1, int p2)
{
    return (pietre[p1].interesse > pietre[p2].interesse);
}


void load_data()
{
    in >> N >> S;                      // 1) N: numero citta' S: posizione della base degli avengers
    in >> M >> C >> R >> Vmin >> Vmax; // 2) M: numero di pietre, C: capacita' del guanto, R:consumo di energia[double], Vmin e Vmax[double]

    if (C != 0)
    {
        Vpond = ((Vmax - Vmin) / (double)C); //controllare che C non sia zero
    }
    else
    {
        Vpond = Vmax;
    }

    grafo.resize(N);
    pietre.resize(M);
    tmpRaccolta.resize(M);
    
    tmpC = C;

    for (int i = 0; i < M; i++) // M) m: massa della pietra i-sima, e: energia della pietra i-sima [int], [int]
    {
        in >> pietre[i].peso >> pietre[i].energia;
        pietre[i].index = i;
        pietre[i].interesse = pietre[i].energia / pietre[i].peso;
        tmpRaccolta[i] = -1;
        pietreInteressanti.push_back(i);
    }

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
            grafo[indexCity].pietreNodo.insert(i);
        }
    }

    int costx;
    for (int i = 1; i < N; i++) // N-1) matrice di adiacenza delle citta, solo parte inferiore con costi
    {
        for (int j = 0; j < i; j++)
        {
            in >> costx;
            grafo[i].vic.push_back(make_pair(j, costx));
            grafo[j].vic.push_back(make_pair(i, costx));
        }
    }
    sort(pietreInteressanti.begin(), pietreInteressanti.end(), sortPietre);
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
    vector<int> pietreSulPercorso;
    pietreSulPercorso.resize(cammino);
    reverse(tmpPercorso.begin(), tmpPercorso.end());

    for (int i = 0; i < cammino; i++) //inizializzazione a -1, nessuna pietra presa
    {
        pietreSulPercorso[i] = -1;
    }
    for (int i = 0; i < M; i++) //salvataggio delle posizioni in cui prendo le pietre, inverso del vettore raccolta
    {
        if (tmpRaccolta[i] != -1)
            pietreSulPercorso[tmpRaccolta[i]] = i;
    }
    for (int i = 0; i < cammino; i++)
    {
        if (pietreSulPercorso[i] != -1)
        {
            PesoAtt += pietre[pietreSulPercorso[i]].peso;
            Eatt += (double)pietre[pietreSulPercorso[i]].energia;
        }
        if (C != 0) // calcolo della' velocita' in uscita dal nodo dopo aver preso la pietra
        {
            Vatt = (double)Vmax - (((double)PesoAtt) * Vpond);
        }
        int vicino;
        for (int j = 0; j < grafo[tmpPercorso[i]].vic.size(); j++) //FIX SCHIFOSO - invece di i era tmpPercorso[i]
        {
            if (grafo[tmpPercorso[i]].vic[j].first == tmpPercorso[i + 1])
                vicino = j;
        }

        Tatt += ((double)(grafo[tmpPercorso[i]].vic[vicino].second)) / Vatt; // T(per un determinato tratto) = lunghezza del tratto/Velocita' per percorrerlo l/v
    }
    G = (double)Eatt;
    T = Tatt;
    E = G - R * T;
    //E(p,t) = G(p) - R*T(p,t)
}

double calcolaCoeffSoloDistanza(int dist, int index)
{
    double coeff = 0.0;
    coeff = ((double)1) / ((double)dist);
    return coeff;
}

void process_solution_greedy()
{
    set<int> nodiToVisit;
    for (int i = 0; i < N; i++) //inserisco tutti i nodi da visitare, tranne la sorgente
    {
        nodiToVisit.insert(i);
    }

    int nodeHere = S;

    while (!nodiToVisit.empty()) //finche ci sono nodi da visitare
    {

        tmpPercorso.push_back(nodeHere);
        grafo[nodeHere].visited = true;
        vector<double> coeff; //vettore dei coefficienti dei vicini del nodo
        coeff.resize(grafo[nodeHere].vic.size());

        if (grafo[nodeHere].pietreNodo.size() > 0 && tmpC > 0)
        {

            int bestPietra=-1;
            int maxRapp = -1;

            for (int i = 0; i < M; i++)
            {
                if (find(grafo[nodeHere].pietreNodo.begin(), grafo[nodeHere].pietreNodo.end(), pietreInteressanti[i]) != grafo[nodeHere].pietreNodo.end())
                {
                    if(pietre[pietreInteressanti[i]].peso <= tmpC)
                    {
                        bestPietra=pietreInteressanti[i];
                        pietreInteressanti.erase(pietreInteressanti.begin() + pietreInteressanti[i]);
                        break;
                    }
                }
            }            
            if(bestPietra != -1) 
            {
                tmpRaccolta[bestPietra] = nodeHere;
                tmpC -= pietre[bestPietra].peso;
            }
        }

        for (int i = 0; i < grafo[nodeHere].vic.size(); i++) //calcola coefficienti dei vicini
        {
            int vicino = grafo[nodeHere].vic[i].first;
            int distanza = grafo[nodeHere].vic[i].second;

            if (vicino == S)
            {
                coeff[i] = INT32_MIN;
                continue;
            }

            if (grafo[vicino].visited)
            {
                coeff[i] = INT32_MIN;
            }
            else
            {
                coeff[i] = calcolaCoeffSoloDistanza(distanza, vicino);
            }
        }

        int massimo = 0;

        for (int i = 0; i < coeff.size(); i++)
        {
            if (coeff[i] > coeff[massimo] && !grafo[grafo[nodeHere].vic[i].first].visited)
            {
                massimo = i;
            }
        }
        massimo = grafo[nodeHere].vic[massimo].first;

        nodiToVisit.erase(nodeHere);
        nodeHere = massimo;
    }
    tmpPercorso.push_back(S);
    get_results();
}

int main()
{
    
    load_data();
    process_solution_greedy();
    print_solution();
    return 0;
}