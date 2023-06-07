#include "empty_class.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <fstream>

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
    Vector2d coordOrigin = mesh.coordinates0D[Vertices1D[0]];
    Vector2d coordEnd= mesh.coordinates0D[Vertices1D[1]];
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
//    swap(vectt[0], vectt[i]);
//    MakeHeap(vectt, i, 0);
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
                idAltroMaxEdge = mesh.vectt[MatriceAdiacenza.Matr[longest][i]].maxedge();
            } 
        }
    }


    // salvo vertici e lati che poi dovrò aggiornare
    array<unsigned int, 3> latiTriAggiornato = triangleToBisect.Edges;
    array<unsigned int, 3> latiTriNuovo = triangleToBisect.Edges;
    array<unsigned int, 3> vertTriAggiornato = triangleToBisect.Vertices2D;
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
        if (vertTriAggiornato[i] != opposite && vertTriAggiornato[i] != mesh.vects[longest].Vertices1D[0]) {
            vertTriAggiornato[i] = newVertex.Id0D;
            break;
        };
    }


    // lati effettivi del triangolo aggiornato
    for (unsigned int i=0; i < 3; i++) {
        if ((mesh.vects[latiTriAggiornato[i]].Vertices1D[0] == opposite && mesh.vects[latiTriAggiornato[i]].Vertices1D[1] == newSegment.Vertices1D[1]) || (mesh.vects[latiTriAggiornato[i]].Vertices1D[1] == opposite && mesh.vects[latiTriAggiornato[i]].Vertices1D[0] == newSegment.Vertices1D[1])) {
            latiTriAggiornato[i] = Mediana.Id1D;
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

    // aggiornamento triangolo vecchio
    triangleToBisect.Edges = latiTriAggiornato;
    triangleToBisect.Vertices2D = vertTriAggiornato;

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

        array<unsigned int, 3> latiTriangoloAggiornato = Triangolo.Edges;
        array<unsigned int, 3> latiUltimoTri = Triangolo.Edges;
        array<unsigned int, 3> vertTriangoloAggiornato = Triangolo.Vertices2D;
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

        // aggiorno
        for (unsigned int i = 0; i<3; i++) {
            if (vertTriangoloAggiornato[i] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) {
                vertTriangoloAggiornato[i] = mesh.vects[idLatoTagliatoNuovo].Vertices1D[0];
                break;
            }
        }

        for (unsigned int i=0; i < 3; i++) {
            if ((mesh.vects[latiTriangoloAggiornato[i]].Vertices1D[0] == newopposite && mesh.vects[latiTriangoloAggiornato[i]].Vertices1D[1] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1]) || (mesh.vects[latiTriangoloAggiornato[i]].Vertices1D[1] == newopposite && mesh.vects[latiTriangoloAggiornato[i]].Vertices1D[0] == mesh.vects[idLatoTagliatoNuovo].Vertices1D[1])) {
                latiTriangoloAggiornato[i] = Unione.Id1D;
            }
        }
        
        Triangolo.Vertices2D = vertTriangoloAggiornato;
        Triangolo.Edges = latiTriangoloAggiornato;
        
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
                break;
            }
        }


    }
    else {
       unsigned int idAltroMaxEdgePropa = NULL;
        unsigned int idAltroTriPropa = NULL;
        
        array<unsigned int, 3> latiTriAggiornatoPropa = Triangolo.Edges;
        array<unsigned int, 3> latiTriNuovoPropa = Triangolo.Edges;
        array<unsigned int, 3> vertTriAggiornatoPropa = Triangolo.Vertices2D;
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
        
        mesh.vects[latoMax].Vertices1D[1] = newVertexPropa.Id0D;
        
        for (unsigned int i = 0;i<3;i++) {
            if (vertTriNuovoPropa[i] != oppositePropa && vertTriNuovoPropa[i] != mesh.vects[latoMax].Vertices1D[1]) {
                vertTriNuovoPropa[i] = newVertexPropa.Id0D;
                break;
            };
        }
    
        for (unsigned int i = 0;i<3;i++) {
            if (vertTriAggiornato[i] != opposite && vertTriAggiornato[i] != mesh.vects[longest].Vertices1D[0]) {
                vertTriAggiornato[i] = newVertex.Id0D;
                break;
            };
        }
        
        for (unsigned int i=0; i < 3; i++) {
            if ((mesh.vects[latiTriAggiornatoPropa[i]].Vertices1D[0] == oppositePropa && mesh.vects[latiTriAggiornatoPropa[i]].Vertices1D[1] == newSegmentPropa.Vertices1D[1]) || (mesh.vects[latiTriAggiornatoPropa[i]].Vertices1D[1] == oppositePropa && mesh.vects[latiTriAggiornatoPropa[i]].Vertices1D[0] == newSegmentPropa.Vertices1D[1])) {
                latiTriAggiornatoPropa[i] = MedianaPropa.Id1D;
                break;
            }
        }
    
        for (unsigned int i=0; i < 3; i++) {
            if (latiTriNuovoPropa[i] == latoMax) {
                latiTriNuovoPropa[i] = newSegmentPropa.Id1D;
            }
            if ((mesh.vects[latiTriNuovoPropa[i]].Vertices1D[0] == oppositePropa && mesh.vects[latiTriNuovoPropa[i]].Vertices1D[1] == mesh.vects[latoMax].Vertices1D[0]) || (mesh.vects[latiTriNuovoPropa[i]].Vertices1D[1] == oppositePropa && mesh.vects[latiTriNuovoPropa[i]].Vertices1D[0] == mesh.vects[latoMax].Vertices1D[0])) {
                latiTriNuovoPropa[i] = MedianaPropa.Id1D;
            }
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        array<unsigned int, 3> latiTriResiduoPropa = Triangolo.Edges;
        array<unsigned int, 3> vertTriResiduoPropa = Triangolo.Vertices2D;
                
        // punto medio vecchio da salvare come origin nella mediana
        // aggiornato = lato tagliato
        // nuovo = altro da prppagazione
        // residuo = residuo da bisezione precedente collegando punti medi
        // taglio lato lungo e poi collego
        // in questo caso è ricorsiva
    }
}


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

mesh.id0D.reserve(mesh.numbercell0D);
mesh.coordinates0D.reserve(mesh.numbercell0D);

for (const string& line : listLines)
{
istringstream converter(line);

unsigned int id;
unsigned int marker;
Vector2d coord;

converter >>  id >> marker >> coord(0) >> coord(1);

mesh.id0D.push_back(id);
mesh.coordinates0D.push_back(coord);
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

  mesh.id1D.reserve(mesh.numbercell1D);
  mesh.vertices1D.reserve(mesh.numbercell1D);

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
    mesh.id1D.push_back(id);
    mesh.vertices1D.push_back(vertices);

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

  mesh.id2D.reserve(mesh.numbercell2D);
  mesh.vertices2D.reserve(mesh.numbercell2D);
  mesh.edges2D.reserve(mesh.numbercell2D);

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
    mesh.id2D.push_back(id);
    mesh.vertices2D.push_back(vertices);
    mesh.edges2D.push_back(edges);
  }
  file.close();
  return true;
}
