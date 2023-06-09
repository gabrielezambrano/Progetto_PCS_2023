#include "empty_class.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <algorithm>

using namespace std;




Cells::TriangularMesh mesh;

bool ImportCell0Ds();
bool ImportCell1Ds();
bool ImportCell2Ds();
void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Cells::Cell2D Triangolo, unsigned int latoMax);

Cells::MatrAdiac MatriceAdiacenza;

namespace Cells
{


// COSTRUTTORI

Cell0D::Cell0D(unsigned int id, unsigned int marker, Vector2d coord)
    {
    unsigned int Id0D = id;
    unsigned int marker0D = marker;
    Vector2d Coord = coord;
    };



Cell1D::Cell1D(unsigned int id,unsigned int marker,Vector2i vertices)
    {
    unsigned int Id1d = id;
    unsigned int marker1D = marker;
    Vector2i Vertices1d = vertices;
    };

Cell2D::Cell2D(unsigned int id,array<unsigned int, 3> Vertices, array<unsigned int, 3> Edges)
    {
    unsigned int Id2D = id;
    array<unsigned int, 3> Vertices2D = Vertices;
    array<unsigned int, 3> Edges2D = Edges;
    };


//METODI LENEDGE, MAXEDGE, (AREA vince)

// !! se nella parte iterativa il lato non viene tolto ma aggiornato con nuova end, non c'è riscalamento id nei vettori -> non serve ciclo for !!

//metterei double anziché void
double Cells::Cell1D::LengthEdge(){
    Vector2d coordOrigin = mesh.vectp[Vertices1D[0]].Coord;
    Vector2d coordEnd= mesh.vectp[Vertices1D[1]].Coord;
    //LengthEdges = (coordEnd-coordOrigin).norm();
    double len = sqrt(pow(coordOrigin[0]-coordEnd[0], 2)+pow(coordOrigin[1] - coordEnd[1], 2));
    return len;
    }


//PROBLEMA TOLLERANZA
unsigned int Cells::Cell2D::maxedge(){
    unsigned int indmax = 0;
    double max = 0.0;
    array<double, 3> lenedges = {0, 0, 0};
    for (unsigned int i = 0; i<3; i++)
    {
        lenedges[i] = mesh.LengthEdges[Edges[i]];
        if(lenedges[i] > max)
        {
            max = lenedges[i];
            indmax = Edges[i]; // oppure indmax = i         NON SO SE HA PIU' SENSO CHE RESTITUISCA L'ID DEL LATO O IL SUO INDICE NEL VETTORE DA TRE
        }
    }
    return indmax;

}



//----------------------------------------------------------------------

// RAFFINAMENTO:

// MATRICE DI ADIACENZA [OK]
// ORDINAMENTO PER AREA [Riccardo]
// PUNTO MEDIO LATO LUNGO (CHE SI AGGIORNA CON STESSA ORIGIN E NUOVA END, LEN = LEN/2) (NEW LATO, LEN = LEN/2)
// TAGLIO LATO LUNGO


MatrAdiac::MatrAdiac() {
    vector<vector<unsigned int>> MatrAdiac(mesh.numbercell1D);
    for (unsigned int i = 0; i < mesh.numbercell2D; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            MatrAdiac[mesh.vectt[i].Edges[j]].push_back(mesh.vectt[i].Id2D);
        }
    }
    Matr = MatrAdiac;
}




//void MakeHeap(vector<Cells::Cell2D> vectt, int n, int i){
//
//  int max = i;
//  int l = 2 * i + 1;
//  int r = 2 * i + 2;
//
//  if (l < n && vectt[l].Area() < vectt[max].Area())
//  {
//    max = l;
//  }
//
//  if (r < n && vectt[r].Area() < vectt[max].Area())
//  {
//    max = r;
//  }
//
//  if (max != i)
//  {
//    swap(vectt[i], vectt[l]);
//
//    MakeHeap(vectt, n, i);
//  }
//}
//
//void HeapSort(vector<Cells::Cell2D>  vectt, int n){
//
//  for (int i = n / 2 - 1; i >= 0; i--)
//    MakeHeap(vectt, n, i);
//
//  for (int i = n - 1; i >= 0; i--)
//  {
//    swap(vectt[0], vectt[i]); // scambio il primo elemento (il più grande) con l'ultimo per escluderlo dal ciclo successivo
//    MakeHeap(vectt, i, 0); // ripeto il ciclo con l'elemento maggiore escluso
//  }
//}



void Bisect(Cells::Cell2D triangleToBisect){

    unsigned int longest = triangleToBisect.maxedge();

    //serve subito controllo: 1) marker lato lungo 2) lato lungo dell'eventuale altro triangolo
    unsigned int markerMaxEdge = mesh.vects[longest].marker1D;
    unsigned int idAltroMaxEdge = NULL;
    unsigned int idAltroTri = NULL;
    if (markerMaxEdge == 0) {
        for (unsigned int i = 0; i<2; i++) {
            if (MatriceAdiacenza.Matr[longest][i] != triangleToBisect.Id2D) {
                idAltroTri = MatriceAdiacenza.Matr[longest][i];
                idAltroMaxEdge = mesh.vectt[idAltroMaxEdge].maxedge();
            }
        }
    }


    // salvo vertici e lati che poi dovrò aggiornare
    array<unsigned int, 3> latiTriNuovo = triangleToBisect.Edges;
    array<unsigned int, 3> vertTriNuovo = triangleToBisect.Vertices2D;
    // inizio bisezione
    Vector2d midCoord;
    midCoord[0] = (mesh.vectp[mesh.vects[longest].Vertices1D[0]].Coord[0] + mesh.vectp[mesh.vects[longest].Vertices1D[1]].Coord[0]) *0.5;
    midCoord[1] = (mesh.vectp[mesh.vects[longest].Vertices1D[0]].Coord[1] + mesh.vectp[mesh.vects[longest].Vertices1D[1]].Coord[1]) *0.5;

    unsigned int markerP;
    if (mesh.vectp[mesh.vects[longest].Vertices1D[0]].marker0D == 0 || mesh.vectp[mesh.vects[longest].Vertices1D[1]].marker0D == 0) {
        markerP = 0;
    }
    else {
    markerP = mesh.vectp[mesh.vects[longest].Vertices1D[0]].marker0D; // per come sono i dati di partenza non ci sono/possono
                                                                      // essere ulteriori configurazioni
    }
    unsigned int newIndexpoint = mesh.vectp.size();
    Cells::Cell0D newVertex = Cell0D(newIndexpoint, markerP, midCoord);
    mesh.vectp.push_back(newVertex);

    unsigned int opposite = NULL;
    for(unsigned int i = 0; i < 3; i++)
    {
        if(!(mesh.vects[longest].Vertices1D[0] == triangleToBisect.Vertices2D[i] || mesh.vects[longest].Vertices1D[1] == triangleToBisect.Vertices2D[i]))
            {
            opposite = triangleToBisect.Vertices2D[i];
            break;
            }
    }

    Vector2i MedianaVert = {opposite, newVertex};

    unsigned int idNewEdge = mesh.vects.size();

    unsigned int markerMediana = 0; // NON PUO' ESSERE ALTRIMENTI



     //Creo segm mediana
    Cell1D Mediana = Cell1D(idNewEdge, markerMediana, MedianaVert);
    mesh.vects.push_back(Mediana);

    Vector2i NewSegVert = {newVertex.Id0D, mesh.vects[longest].Vertices1D[1]};



    //Creo segm pto medio -> end vecchia
    Cell1D newSegment = Cell1D(idNewEdge + 1, mesh.vects[longest].marker1D, NewSegVert);
    mesh.vects.push_back(newSegment);


    //aggiorno segm origin vecchia -> pto medio
    mesh.vects[longest].Vertices1D[1] = newVertex.Id0D;  // GUARDATO FINO A QUA



    //vertici effettivi del triangolo nuovo
    for (unsigned int i = 0;i<3;i++) {
        if (vertTriNuovo[i] != opposite && vertTriNuovo[i] != mesh.vects[longest].Vertices1D[1]) {
            vertTriNuovo[i] = newVertex.Id0D;
            break;
        };
    }

    // vertici effettivi del triangolo aggiornato
    for (unsigned int i = 0;i<3;i++) {
        if (triangleToBisect.Vertices2D[i] != opposite && triangleToBisect.Vertices2D[i] != mesh.vects[longest].Vertices1D[0]) {
            triangleToBisect.Vertices2D[i] = newVertex.Id0D;
            break;
        };
    }


    // lati effettivi del triangolo aggiornato
    for (unsigned int i=0; i < 3; i++) {
        if ((mesh.vects[triangleToBisect.Edges[i]].Vertices1D[0] == opposite && mesh.vects[triangleToBisect.Edges[i]].Vertices1D[1] == newSegment.Vertices1D[1]) || (mesh.vects[triangleToBisect.Edges[i]].Vertices1D[1] == opposite && mesh.vects[triangleToBisect.Edges[i]].Vertices1D[0] == newSegment.Vertices1D[1])) {
            triangleToBisect.Edges[i] = Mediana.Id1D;
            break;
        }
    }

    // lati effettivi del triangolo nuovo


    for (unsigned int i=0; i < 3; i++) {
        if (latiTriNuovo[i] == longest) {
            latiTriNuovo[i] = newSegment.Id1D;
        }
        if ((mesh.vects[latiTriNuovo[i]].Vertices1D[0] == opposite && mesh.vects[latiTriNuovo[i]].Vertices1D[1] == mesh.vects[longest].Vertices1D[0]) || (mesh.vects[latiTriNuovo[i]].Vertices1D[1] == opposite && mesh.vects[latiTriNuovo[i]].Vertices1D[0] == mesh.vects[longest].Vertices1D[0])) {
            latiTriNuovo[i] = Mediana.Id1D;
        }
    }

    // creazione triangolo nuovo
    unsigned int idnewTri = mesh.vectt.size();
    Cell2D newTri = Cell2D(idnewTri, vertTriNuovo, latiTriNuovo);
    mesh.vectt.push_back(newTri);


    // aggiorno matrice di adiacenza
    // aggiungo mediana
    MatriceAdiacenza.Matr.push_back({newTri.Id2D, triangleToBisect.Id2D});
    // aggiungo newSegment
    MatriceAdiacenza.Matr.push_back({newTri.Id2D, idAltroTri});
    // aggiorno terzo lato
    for (unsigned int i=0; i<3; i++) {
        if(newTri.Edges[i] != Mediana.Id1D && newTri.Edges[i] != newSegment.Id1D){
            for (unsigned int j = 0; j < 2; j++) {
                if(MatriceAdiacenza.Matr[newTri.Edges[i]][j] == triangleToBisect.Id2D){
                    MatriceAdiacenza.Matr[newTri.Edges[i]][j] = newTri.Id2D;
                    break;
                }
            }
            break;
        }
    }

    if (markerMaxEdge == 0) {
        Propagazione(longest, newSegment.Id1D, mesh.vectt[idAltroTri], idAltroMaxEdge);
    }


} // fine Bisect



void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Cell2D Triangolo, unsigned int latoMax){
    if (idLatoTagliatoVecchio == latoMax){

        // collega pto medio e vertice opposto

        array<unsigned int, 3> latiUltimoTri = Triangolo.Edges;
        array<unsigned int, 3> vertUltimoTri = Triangolo.Vertices2D;

        unsigned int newopposite = NULL;
        for(unsigned int i = 0; i < 3; i++)
        {
            if(!(mesh.vects[idLatoTagliatoVecchio].Vertices1D[0] == Triangolo.Vertices2D[i] || mesh.vects[idLatoTagliatoNuovo].Vertices1D[1] == Triangolo.Vertices2D[i]))
                {
                newopposite = Triangolo.Vertices2D[i];
                break;
                }
        }
        //creo ultimo lato
        Cell1D Unione = Cell1D(mesh.vectp.size(), 0, {mesh.vects[idLatoTagliatoNuovo].Vertices1D[0], newopposite});
        mesh.vects.push_back(Unione);

        // creo Ultimo triangolo
        for (unsigned int i = 0; i<3; i++) {
            if (vertUltimoTri[i] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                vertUltimoTri[i] = mesh.vects[idLatoTagliatoVecchio].Vertices1D[1];
                break;
            }
        }

        for (unsigned int i=0; i < 3; i++) {
            if (latiUltimoTri[i] == latoMax) {
                latiUltimoTri[i] = idLatoTagliatoNuovo;
            }
            if ((mesh.vects[latiUltimoTri[i]].Vertices1D[0] == newopposite && mesh.vects[latiUltimoTri[i]].Vertices1D[1] == mesh.vects[latoMax].Vertices1D[0]) || (mesh.vects[latiUltimoTri[i]].Vertices1D[1] == newopposite && mesh.vects[latiUltimoTri[i]].Vertices1D[0] == mesh.vects[latoMax].Vertices1D[0])) {
                latiUltimoTri[i] = Unione.Id1D;
            }
        }

        Cell2D UltimoTriangolo = Cell2D(mesh.vects.size(), vertUltimoTri, latiUltimoTri);
        mesh.vectt.push_back(UltimoTriangolo);

        // aggiorno vertici
        for (unsigned int i = 0; i<3; i++) {
            if (Triangolo.Vertices2D[i] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
                Triangolo.Vertices2D[i] = mesh.vects[idLatoTagliatoNuovo].Vertices1D[0];
                break;
            }
        }

        for (unsigned int i=0; i < 3; i++) {
            if ((mesh.vects[Triangolo.Edges[i]].Vertices1D[0] == newopposite && mesh.vects[Triangolo.Edges[i]].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) || (mesh.vects[Triangolo.Edges[i]].Vertices1D[1] == newopposite && mesh.vects[Triangolo.Edges[i]].Vertices1D[0] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1])) {
                Triangolo.Edges[i] = Unione.Id1D;
            }
        }


        //aggiorno matrice di adiacenza
        MatriceAdiacenza.Matr.push_back({Triangolo.Id2D, UltimoTriangolo.Id2D});


        for (unsigned int i=0; i<3; i++) {
            if(UltimoTriangolo.Edges[i] != Unione.Id1D && UltimoTriangolo.Edges[i] != idLatoTagliatoNuovo){
                for (unsigned int j = 0; j < 2; j++) {
                    if(MatriceAdiacenza.Matr[UltimoTriangolo.Edges[i]][j] == Triangolo.Id2D){
                        MatriceAdiacenza.Matr[UltimoTriangolo.Edges[i]][j] = UltimoTriangolo.Id2D;
                        break;
                    }
                }
                
            }
                        
        }
        for (unsigned int i=0; i<2; i++) {
            if (MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] = UltimoTriangolo.Id2D;
            }
        }


    }
    else {


        // va fatto un if che mette l'end di latoMax = end latotagliatonuovo (caso 1)
        //                            origin latoMax = origin latotagliatovecchio (caso 2)

        for (unsigned int i = 0; i <2 ; i++) {
            if (mesh.vects[latoMax].Vertices1D[i] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0] && 1==1){
                unsigned int supporto = mesh.vects[latoMax].Vertices1D[i];
                mesh.vects[latoMax].Vertices1D[i] = mesh.vects[latoMax].Vertices1D[i-1];
                mesh.vects[latoMax].Vertices1D[i-1] = supporto;

            }
            else if (mesh.vects[latoMax].Vertices1D[i] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[1] && 1==0) {
                unsigned int supporto = mesh.vects[latoMax].Vertices1D[i];
                mesh.vects[latoMax].Vertices1D[i] = mesh.vects[latoMax].Vertices1D[i+1];
                mesh.vects[latoMax].Vertices1D[i+1] = supporto;
            }
        }



        unsigned int markerMaxEdgePropa = mesh.vects[latoMax].marker1D;
        unsigned int idAltroMaxEdgePropa = NULL;
        unsigned int idAltroTriPropa = NULL;
        if (markerMaxEdgePropa == 0) {
            for (unsigned int i = 0; i<2; i++) {
                if (MatriceAdiacenza.Matr[latoMax][i] != Triangolo.Id2D) {
                    idAltroTriPropa = MatriceAdiacenza.Matr[latoMax][i];
                    idAltroMaxEdgePropa = mesh.vectt[MatriceAdiacenza.Matr[latoMax][i]].maxedge();
                }
            }
        }



        array<unsigned int, 3> latiTriNuovoPropa = Triangolo.Edges;
        array<unsigned int, 3> vertTriNuovoPropa = Triangolo.Vertices2D;

        Vector2d midCoordPropa;
        midCoordPropa[0] = (mesh.vectp[mesh.vects[latoMax].Vertices1D[0]].Coord[0] + mesh.vectp[mesh.vects[latoMax].Vertices1D[1]].Coord[0]) *0.5;
        midCoordPropa[1] = (mesh.vectp[mesh.vects[latoMax].Vertices1D[0]].Coord[1] + mesh.vectp[mesh.vects[latoMax].Vertices1D[1]].Coord[1]) *0.5;

        unsigned int markerPPropa;
        if (mesh.vectp[mesh.vects[latoMax].Vertices1D[0]].marker0D == 0 || mesh.vectp[mesh.vects[latoMax].Vertices1D[1]].marker0D == 0) {
            markerPPropa = 0;
        }
        else {
        markerPPropa = mesh.vectp[mesh.vects[latoMax].Vertices1D[0]].marker0D; // per come sono i dati di partenza non ci sono/possono
                                                                               // essere ulteriori configurazioni
        }

        unsigned int newIndexpointPropa = mesh.vectp.size();
        Cells::Cell0D newVertexPropa = Cell0D(newIndexpointPropa, markerPPropa, midCoordPropa);
        mesh.vectp.push_back(newVertexPropa);

        unsigned int oppositePropa = NULL;
        for(unsigned int i = 0; i < 3; i++)
        {
            if(!(mesh.vects[latoMax].Vertices1D[0] == Triangolo.Vertices2D[i] || mesh.vects[latoMax].Vertices1D[1] == Triangolo.Vertices2D[i]))
                {
                oppositePropa = Triangolo.Vertices2D[i];
                break;
                }
        }

        Vector2i MedianaVertPropa = {oppositePropa, newVertexPropa};

        unsigned int idNewEdgePropa = mesh.vects.size();

        unsigned int markerMedianaPropa = 0; // NON PUO' ESSERE ALTRIMENTI

        Cell1D MedianaPropa = Cell1D(idNewEdgePropa, markerMedianaPropa, MedianaVertPropa);
        mesh.vects.push_back(MedianaPropa);

        Vector2i NewSegVertPropa = {newVertexPropa.Id0D, mesh.vects[latoMax].Vertices1D[1]};


        Cell1D newSegmentPropa = Cell1D(idNewEdgePropa + 1, mesh.vects[latoMax].marker1D, NewSegVertPropa);
        mesh.vects.push_back(newSegmentPropa);

        // aggiorno end pto medio propa
        mesh.vects[latoMax].Vertices1D[1] = newVertexPropa.Id0D; //GUARDATO FINO A QUI


        unsigned int lontano = NULL;
        for (unsigned int i = 0;i<3;i++) {
            if (Triangolo.Vertices2D[i] != mesh.vects[idLatoTagliatoNuovo].Vertices1D[1] && Triangolo.Vertices2D[i] != mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                lontano = Triangolo.Vertices2D[i];
            }
        }

        // vertici effettivi nuovo propa
        for (unsigned int i = 0;i<3;i++) {
            if (vertTriNuovoPropa[i] == lontano) {
                vertTriNuovoPropa[i] = newVertexPropa.Id0D;
                break;
            }
        } // ok

        // vertici aggiornati
        for (unsigned int i = 0;i<3;i++) {
            if (Triangolo.Vertices2D[i] != oppositePropa && Triangolo.Vertices2D[i] != lontano) {
                Triangolo.Vertices2D[i] = newVertexPropa.Id0D;
                break;
            };
        } // ok

        // lati aggiornati
        if (mesh.vects[latoMax].Vertices1D[0] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            for (unsigned int i=0; i < 3; i++) {
                if (Triangolo.Edges[i] == latoMax){
                    Triangolo.Edges[i] = newSegmentPropa.Vertices1D[1];
                }
                if (Triangolo.Edges[i] == idLatoTagliatoVecchio){
                    Triangolo.Edges[i] = MedianaPropa.Id1D;
                }
            }

        }
        else if (mesh.vects[latoMax].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i=0; i < 3; i++) {
                if (Triangolo.Edges[i] == idLatoTagliatoVecchio) {
                    Triangolo.Edges[i] = MedianaPropa.Id1D;
                }
            }
        }

        // lati effettivi nuovo propa
        if (mesh.vects[latoMax].Vertices1D[0] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            for (unsigned int i=0; i < 3; i++) {
                if (latiTriNuovoPropa[i] != latoMax && latiTriNuovoPropa[i] != idLatoTagliatoVecchio){
                    latiTriNuovoPropa[i] = MedianaPropa.Id1D;
                }
            }
        }
        else if (mesh.vects[latoMax].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i=0; i < 3; i++) {
                if (latiTriNuovoPropa[i] == latoMax){
                    latiTriNuovoPropa[i] = newSegmentPropa.Id1D;
                }
                if ((mesh.vects[latiTriNuovoPropa[i]].Vertices1D[0] == oppositePropa && mesh.vects[latiTriNuovoPropa[i]].Vertices1D[1] == lontano) || (mesh.vects[latiTriNuovoPropa[i]].Vertices1D[1] == oppositePropa && mesh.vects[latiTriNuovoPropa[i]].Vertices1D[0] == lontano)) {
                    latiTriNuovoPropa[i] = MedianaPropa.Id1D;
                }
            }
        }

        Cell2D Penultimo = Cell2D(mesh.vectt.size(), vertTriNuovoPropa, latiTriNuovoPropa);
        mesh.vectt.push_back(Penultimo);
        
        //aggiorno PARZIALMENTE matrice di adiac
        
        if (mesh.vects[latoMax].Vertices1D[0] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            MatriceAdiacenza.Matr.push_back({Triangolo.Id2D, Penultimo.Id2D});
            MatriceAdiacenza.Matr.push_back({Triangolo.Id2D, idAltroTriPropa});
            for (unsigned int i = 0; i<2; i++) {
                if (MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                    MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] = Penultimo.Id2D;
                }
            }
        }
        else if (mesh.vects[latoMax].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            MatriceAdiacenza.Matr.push_back({Triangolo.Id2D, Penultimo.Id2D});
            for (unsigned int i = 0; i<2; i++) {
                if (MatriceAdiacenza.Matr[idLatoTagliatoVecchio][i] == Triangolo.Id2D) {
                    MatriceAdiacenza.Matr[idLatoTagliatoVecchio][i] = Penultimo.Id2D;
                }
                
            }
        }
        
        array<unsigned int, 3> latiTriResiduoPropa = Penultimo.Edges;
        array<unsigned int, 3> vertTriResiduoPropa = Penultimo.Vertices2D;
        
        
        Cell1D Unione = Cell1D(mesh.vects.size(), 0, {mesh.vects[idLatoTagliatoVecchio].Vertices1D[1], newVertexPropa});
        mesh.vects.push_back(Unione);
        
        if (mesh.vects[latoMax].Vertices1D[0] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            
            for (unsigned int i = 0; i<3; i++) {
                // aggiorno vertici
                if (Penultimo.Vertices2D[i] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                    Penultimo.Vertices2D[i] = mesh.vects[idLatoTagliatoVecchio].Vertices1D[1];
                }
                if (vertTriResiduoPropa[i] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
                    vertTriResiduoPropa[i] = mesh.vects[idLatoTagliatoNuovo].Vertices1D[0];
                }
                // aggiorno lati
                if (Penultimo.Edges[i] == idLatoTagliatoVecchio) {
                    Penultimo.Edges[i] = idLatoTagliatoNuovo;
                }
                if (Penultimo.Edges[i] == latoMax) {
                    Penultimo.Edges[i] = Unione.Id1D;
                }
                if (latiTriResiduoPropa[i] == MedianaPropa.Id1D) {
                    latiTriResiduoPropa[i] = Unione.Id1D;
                }
                
            }
            Cell2D Residuo = Cell2D(mesh.vectt.size(), vertTriResiduoPropa, latiTriResiduoPropa);
            mesh.vectt.push_back(Residuo);
            
            // aggiorno matrice adiacenza
            MatriceAdiacenza.Matr.push_back({Penultimo.Id2D, Residuo.Id2D});
            for (unsigned int i = 0; i<2; i++) {
                if (MatriceAdiacenza.Matr[latoMax][i] == Triangolo.Id2D) {
                    MatriceAdiacenza.Matr[latoMax][i] = Residuo.Id2D;
                }
                if (MatriceAdiacenza.Matr[idLatoTagliatoVecchio][i] == Triangolo.Id2D) {
                    MatriceAdiacenza.Matr[idLatoTagliatoVecchio][i] = Residuo.Id2D;
                }
            }
        }
        else if (mesh.vects[latoMax].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i = 0; i<3; i++) {
                if (Penultimo.Vertices2D[i] == newSegmentPropa.Vertices1D[1]) {
                    Penultimo.Vertices2D[i] = mesh.vects[idLatoTagliatoVecchio].Vertices1D[1];
                }
                if (vertTriResiduoPropa[i] == mesh.vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                    vertTriResiduoPropa[i] = mesh.vects[idLatoTagliatoNuovo].Vertices1D[0];
                }
                // aggiorno lati
                if (Penultimo.Edges[i] == newSegmentPropa.Id1D) {
                    Penultimo.Edges[i] = Unione.Id1D;
                }
                if (latiTriResiduoPropa[i] == MedianaPropa.Id1D) {
                    latiTriResiduoPropa[i] = Unione.Id1D;
                }
                if (latiTriResiduoPropa[i] == idLatoTagliatoVecchio) {
                    latiTriResiduoPropa[i] = idLatoTagliatoNuovo;
                }
            }
            
            Cell2D Residuo = Cell2D(mesh.vectt.size(), vertTriResiduoPropa, latiTriResiduoPropa);
            mesh.vectt.push_back(Residuo);
            
            // aggiorno matrice adiacenza
            MatriceAdiacenza.Matr.push_back({Residuo.Id2D, idAltroTriPropa});
            MatriceAdiacenza.Matr.push_back({Penultimo.Id2D, Residuo.Id2D});
            for (unsigned int i = 0; i<2; i++) {                
                if (MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                    MatriceAdiacenza.Matr[idLatoTagliatoNuovo][i] = Residuo.Id2D;
                }
            }
        }
        
        // void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Cell2D Triangolo, unsigned int latoMax){
        // in questo caso è ricorsiva
        if (markerMaxEdgePropa == 0) {
            Propagazione(latoMax, newSegmentPropa.Id1D, mesh.vectt[idAltroTriPropa], idAltroMaxEdgePropa);
        }
    } // fine else (lato lungo diverso dal precedente)
} // fine Propagazione


} // fine namespace Cells


//-----------------------------------------------------------------------

//IMPORTAZIONE

bool ImportCell0Ds()
{

ifstream file;
file.open("./Cell0Ds.csv");

if(file.fail())
return false;

list<string> listLines;
string line;
while (getline(file, line))
listLines.push_back(line);

file.close();

listLines.pop_front();

mesh.numbercell0D = listLines.size();

if (mesh.numbercell0D == 0)
{
cerr << "There is no cell 0D" << endl;
return false;
}



for (const string& line : listLines)
{
istringstream converter(line);

unsigned int id;
unsigned int marker;
Vector2d coord;

converter >>  id >> marker >> coord(0) >> coord(1);


Cells::Cell0D point = Cells::Cell0D(id,marker,coord);
mesh.vectp.push_back(point);

//    if( marker != 0)
//    {
//      if (mesh.Cell0DMarkers.find(marker) == mesh.Cell0DMarkers.end())
//        mesh.Cell0DMarkers.insert({marker, {id}});
//      else
//        mesh.Cell0DMarkers[marker].push_back(id);
//    }

  }
  file.close();
  return true;
}

bool ImportCell1Ds()
{

  ifstream file;
  file.open("./Cell1Ds.csv");

  if(file.fail())
    return false;

  list<string> listLines;
  string line;
  while (getline(file, line))
    listLines.push_back(line);

  listLines.pop_front();

  mesh.numbercell1D = listLines.size();

  if (mesh.numbercell1D == 0)
  {
    cerr << "There is no cell 1D" << endl;
    return false;
  }



  for (const string& line : listLines)
  {
    istringstream converter(line);

    unsigned int id;
    unsigned int marker;
    Vector2i vertices;

    converter >>  id >> marker >> vertices(0) >> vertices(1);
    Cells::Cell1D segment = Cells::Cell1D(id,marker,vertices);
    mesh.vects.push_back(segment);
    mesh.LengthEdges.push_back(segment.LengthEdge());


//    if( marker != 0)
//    {
//      if (mesh.Cell1DMarkers.find(marker) == mesh.Cell1DMarkers.end())
//        mesh.Cell1DMarkers.insert({marker, {id}});
//      else
//        mesh.Cell1DMarkers[marker].push_back(id);
//    }
  }

  file.close();

  return true;
}



bool ImportCell2Ds()
{

  ifstream file;
  file.open("./Cell2Ds.csv");

  if(file.fail())
    return false;

  list<string> listLines;
  string line;
  while (getline(file, line))
    listLines.push_back(line);

  listLines.pop_front();


  mesh.numbercell2D = listLines.size();

  if (mesh.numbercell2D == 0)
  {
    cerr << "There is no cell 2D" << endl;
    return false;
  }



  for (const string& line : listLines)
  {
    istringstream converter(line);

    unsigned int id;
    array<unsigned int, 3> vertices;
    array<unsigned int, 3> edges;

    converter >>  id;
    for(unsigned int i = 0; i < 3; i++)
      converter >> vertices[i];
    for(unsigned int i = 0; i < 3; i++)
      converter >> edges[i];

    Cells::Cell2D triangle = Cells::Cell2D(id,vertices,edges);
    mesh.vectt.push_back(triangle);

  }
  file.close();
  return true;
}
