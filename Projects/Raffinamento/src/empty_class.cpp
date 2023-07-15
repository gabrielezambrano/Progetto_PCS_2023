#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <algorithm>
#include "cmath"
#include "empty_class.hpp"
#include "sorting.hpp"


using namespace std;
using namespace Eigen;
using namespace Sorting;
using namespace Project;




Project::TriangularMesh mesh = Project::TriangularMesh();


//void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Project::Cell2D Triangolo, unsigned int latoMax);


Cell0D::Cell0D(unsigned int id, unsigned int marker, Vector2d coord)
    {
    unsigned int Id0D = id;
    unsigned int marker0D = marker;
    Vector2d Coord = coord;
    };



Cell1D::Cell1D(unsigned int id,unsigned int marker,vector<unsigned int>& vertices)
    {
    unsigned int Id1d = id;
    unsigned int marker1D = marker;
    vector<unsigned int> Vertices1d = vertices;
    };

Cell2D::Cell2D(unsigned int& id,array<unsigned int, 3>& Vertices, array<unsigned int, 3>& Edges)
    {
    unsigned int Id2D = id;
    array<unsigned int, 3> Vertices2D = Vertices;
    array<unsigned int, 3> Edges2D = Edges;
    };







double Project::Cell1D::LengthEdge(vector<Project::Cell0D>& vectp){
    Vector2d coordOrigin = mesh.vectp[this->Vertices1D[0]].Coord;
    Vector2d coordEnd= mesh.vectp[this->Vertices1D[1]].Coord;
    //double len = (coordEnd-coordOrigin).norm();
    double len = sqrt(pow(coordOrigin.x() - coordEnd.x(), 2)+pow(coordOrigin.y() - coordEnd.y(), 2));
    return len;
    }

unsigned int Project::Cell2D::maxedge(vector<Project::Cell0D>& vectp,vector<Project::Cell1D>& vects){ 
    unsigned int indmax = 0;
    double max = mesh.vects[this->Edges[0]].LengthEdge();
    for (unsigned int i = 1; i<3; i++){
        if(mesh.vects[this->Edges[i]].LengthEdge()  > max - tol1D){ 
            max = mesh.vects[this->Edges[i]].LengthEdge();
            indmax = i;
        }
    }
    return this->Edges[indmax];

}

double Project::Cell2D::Area(vector<Project::Cell0D>& vectp){
             //Formula dell'area di Gauss
             double A_12 = (mesh.vectp[this->Vertices2D[0]].Coord[0]*mesh.vectp[this->Vertices2D[1]].Coord[1]) - (mesh.vectp[this->Vertices2D[0]].Coord[1]*mesh.vectp[this->Vertices2D[1]].Coord[0]);
             double A_23 = (mesh.vectp[this->Vertices2D[1]].Coord[0]*mesh.vectp[this->Vertices2D[2]].Coord[1]) - (mesh.vectp[this->Vertices2D[1]].Coord[1]*mesh.vectp[this->Vertices2D[2]].Coord[0]);
             double A_31 = (mesh.vectp[this->Vertices2D[2]].Coord[0]*mesh.vectp[this->Vertices2D[0]].Coord[1]) - (mesh.vectp[this->Vertices2D[2]].Coord[1]*mesh.vectp[this->Vertices2D[0]].Coord[0]);
             return abs((A_12+A_23+A_31)/2);
  }







// !! se nella parte iterativa il lato non viene tolto ma aggiornato con nuova end, non c'è riscalamento id nei vettori -> non serve ciclo for !!

//metterei double anziché void




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


Project::Cell0D point = Project::Cell0D(id,marker,coord);
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
    vector<unsigned int> vertices;

    converter >>  id >> marker >> vertices[0] >> vertices[1];
    Project::Cell1D segment = Project::Cell1D(id,marker,vertices);
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

    Project::Cell2D triangle = Project::Cell2D(id,vertices,edges);
    mesh.vectt.push_back(triangle);

  }
  file.close();
  return true;
}





//----------------------------------------------------------------------

// RAFFINAMENTO:

// MATRICE DI ADIACENZA [OK]
// ORDINAMENTO PER AREA 
// PUNTO MEDIO LATO LUNGO (CHE SI AGGIORNA CON STESSA ORIGIN E NUOVA END, LEN = LEN/2) (NEW LATO, LEN = LEN/2)
// TAGLIO LATO LUNGO


MatrAdiac::MatrAdiac(vector<Project::Cell2D> vectt, vector<Project::Cell1D> vects) {
    vector<vector<unsigned int>> MatrAdiac(vects.size(), vector<unsigned int>());
    for (unsigned int i = 0; i < vectt.size(); i++) {
        for (unsigned int j = 0; j < 3; j++) {
            MatrAdiac[vectt[i].Edges[j]].push_back(vectt[i].Id2D);
        }
    }
    Matr = MatrAdiac;
}





void Bisect(Project::Cell2D triangleToBisect, vector<Project::Cell0D> vectp, vector<Project::Cell1D> vects, vector<Project::Cell2D> vectt, vector<vector<unsigned int>> Matr){


    unsigned int longest = triangleToBisect.maxedge();

    //serve subito controllo: 1) marker lato lungo 2) lato lungo dell'eventuale altro triangolo
    unsigned int markerMaxEdge = vects[longest].marker1D;
    unsigned int idAltroMaxEdge = 0;
    unsigned int idAltroTri = 0;
    if (markerMaxEdge == 0) {
        for (unsigned int i = 0; i<2; i++) {
            if (Matr[longest][i] != triangleToBisect.Id2D) {
                idAltroTri = Matr[longest][i];
                idAltroMaxEdge = vectt[idAltroMaxEdge].maxedge();
            }
        }
    }


    // salvo vertici e lati che poi dovrò aggiornare
    array<unsigned int, 3> latiTriNuovo = triangleToBisect.Edges;
    array<unsigned int, 3> vertTriNuovo = triangleToBisect.Vertices2D;
    // inizio bisezione
    Vector2d midCoord;
    midCoord[0] = (vectp[vects[longest].Vertices1D[0]].Coord[0] + vectp[vects[longest].Vertices1D[1]].Coord[0]) *0.5;
    midCoord[1] = (vectp[vects[longest].Vertices1D[0]].Coord[1] + vectp[vects[longest].Vertices1D[1]].Coord[1]) *0.5;

    unsigned int markerP;
    if (vectp[vects[longest].Vertices1D[0]].marker0D == 0 || vectp[vects[longest].Vertices1D[1]].marker0D == 0) {
        markerP = 0;
    }
    else {
    markerP = vectp[vects[longest].Vertices1D[0]].marker0D; // per come sono i dati di partenza non ci sono/possono
                                                                      // essere ulteriori configurazioni
    }
    unsigned int newIndexpoint = vectp.size();
    Project::Cell0D newVertex = Cell0D(newIndexpoint, markerP, midCoord);
    vectp.push_back(newVertex);

    unsigned int opposite = 0;
    for(unsigned int i = 0; i < 3; i++)
    {
        if(!(vects[longest].Vertices1D[0] == triangleToBisect.Vertices2D[i] || vects[longest].Vertices1D[1] == triangleToBisect.Vertices2D[i]))
            {
            opposite = triangleToBisect.Vertices2D[i];
            break;
            }
    }

    vector<unsigned int> MedianaVert = {opposite, newVertex.Id0D};

    unsigned int idNewEdge = vects.size();
    unsigned int idNewEdge2 = idNewEdge + 1;

    unsigned int markerMediana = 0; // NON PUO' ESSERE ALTRIMENTI



     //Creo segm mediana
    Cell1D Mediana = Cell1D(idNewEdge, markerMediana, MedianaVert);
    vects.push_back(Mediana);

    vector<unsigned int> NewSegVert = {newVertex.Id0D, vects[longest].Vertices1D[1]};



    //Creo segm pto medio -> end vecchia 
    Cell1D newSegment = Cell1D(idNewEdge2, vects[longest].marker1D, NewSegVert);
    vects.push_back(newSegment);


    //aggiorno segm origin vecchia -> pto medio
    vects[longest].Vertices1D[1] = newVertex.Id0D;  // GUARDATO FINO A QUA



    //vertici effettivi del triangolo nuovo
    for (unsigned int i = 0;i<3;i++) {
        if (vertTriNuovo[i] != opposite && vertTriNuovo[i] != vects[longest].Vertices1D[1]) {
            vertTriNuovo[i] = newVertex.Id0D;
            break;
        };
    }

    // vertici effettivi del triangolo aggiornato
    for (unsigned int i = 0;i<3;i++) {
        if (triangleToBisect.Vertices2D[i] != opposite && triangleToBisect.Vertices2D[i] != vects[longest].Vertices1D[0]) {
            triangleToBisect.Vertices2D[i] = newVertex.Id0D;
            break;
        };
    }


    // lati effettivi del triangolo aggiornato
    for (unsigned int i=0; i < 3; i++) {
        if ((vects[triangleToBisect.Edges[i]].Vertices1D[0] == opposite && vects[triangleToBisect.Edges[i]].Vertices1D[1] == newSegment.Vertices1D[1]) || (vects[triangleToBisect.Edges[i]].Vertices1D[1] == opposite && vects[triangleToBisect.Edges[i]].Vertices1D[0] == newSegment.Vertices1D[1])) {
            triangleToBisect.Edges[i] = Mediana.Id1D;
            break;
        }
    }

    // lati effettivi del triangolo nuovo


    for (unsigned int i=0; i < 3; i++) {
        if (latiTriNuovo[i] == longest) {
            latiTriNuovo[i] = newSegment.Id1D;
        }
        if ((vects[latiTriNuovo[i]].Vertices1D[0] == opposite && vects[latiTriNuovo[i]].Vertices1D[1] == vects[longest].Vertices1D[0]) || (vects[latiTriNuovo[i]].Vertices1D[1] == opposite && vects[latiTriNuovo[i]].Vertices1D[0] == vects[longest].Vertices1D[0])) {
            latiTriNuovo[i] = Mediana.Id1D;
        }
    }

    // creazione triangolo nuovo
    unsigned int idnewTri = vectt.size();
    Cell2D newTri = Cell2D(idnewTri, vertTriNuovo, latiTriNuovo);
    vectt.push_back(newTri);


    // aggiorno matrice di adiacenza
    // aggiungo mediana
    Matr.push_back({newTri.Id2D, triangleToBisect.Id2D});
    // aggiungo newSegment
    Matr.push_back({newTri.Id2D, idAltroTri});
    // aggiorno terzo lato
    for (unsigned int i=0; i<3; i++) {
        if(newTri.Edges[i] != Mediana.Id1D && newTri.Edges[i] != newSegment.Id1D){
            for (unsigned int j = 0; j < 2; j++) {
                if(Matr[newTri.Edges[i]][j] == triangleToBisect.Id2D){
                    Matr[newTri.Edges[i]][j] = newTri.Id2D;
                    break;
                }
            }
            break;
        }
    }

    if (markerMaxEdge == 0) {
        Cell2D AltroTri = vectt[idAltroTri];
        Propagazione(longest, newSegment.Id1D, AltroTri, idAltroMaxEdge, vectp, vects, vectt, Matr);
    }


} // fine Bisect



void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Cell2D Triangolo, unsigned int latoMax, vector<Project::Cell0D> vectp, vector<Project::Cell1D> vects, vector<Project::Cell2D> vectt, vector<vector<unsigned int>> Matr){

    if (idLatoTagliatoVecchio == latoMax){ //Parte Umberto

        // collega pto medio e vertice opposto

        array<unsigned int, 3> latiUltimoTri = Triangolo.Edges;
        array<unsigned int, 3> vertUltimoTri = Triangolo.Vertices2D;

        unsigned int newopposite = 0;
        for(unsigned int i = 0; i < 3; i++)
        {
            if(!(vects[idLatoTagliatoVecchio].Vertices1D[0] == Triangolo.Vertices2D[i] || vects[idLatoTagliatoNuovo].Vertices1D[1] == Triangolo.Vertices2D[i]))
                {
                newopposite = Triangolo.Vertices2D[i];
                break;
                }
        }

        unsigned int marker0 = 0;
        unsigned int size = vects.size();
        vector<unsigned int> vettore = {vects[idLatoTagliatoNuovo].Vertices1D[0], newopposite};
        //creo ultimo lato
        Cell1D Unione = Cell1D(size, marker0, vettore);
        vects.push_back(Unione);

        // creo Ultimo triangolo
        for (unsigned int i = 0; i<3; i++) {
            if (vertUltimoTri[i] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                vertUltimoTri[i] = vects[idLatoTagliatoVecchio].Vertices1D[1];
                break;
            }
        }

        for (unsigned int i=0; i < 3; i++) {
            if (latiUltimoTri[i] == latoMax) {
                latiUltimoTri[i] = idLatoTagliatoNuovo;
            }
            if ((vects[latiUltimoTri[i]].Vertices1D[0] == newopposite && vects[latiUltimoTri[i]].Vertices1D[1] == vects[latoMax].Vertices1D[0]) || (vects[latiUltimoTri[i]].Vertices1D[1] == newopposite && vects[latiUltimoTri[i]].Vertices1D[0] == vects[latoMax].Vertices1D[0])) {
                latiUltimoTri[i] = Unione.Id1D;
            }
        }

        unsigned int sizeT = vectt.size();

        Cell2D UltimoTriangolo = Cell2D(sizeT, vertUltimoTri, latiUltimoTri);
        vectt.push_back(UltimoTriangolo);

        // aggiorno vertici
        for (unsigned int i = 0; i<3; i++) {
            if (Triangolo.Vertices2D[i] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
                Triangolo.Vertices2D[i] = vects[idLatoTagliatoNuovo].Vertices1D[0];
                break;
            }
        }

        // aggiorno lati
        for (unsigned int i=0; i < 3; i++) {
            if ((vects[Triangolo.Edges[i]].Vertices1D[0] == newopposite && vects[Triangolo.Edges[i]].Vertices1D[1] == vects[idLatoTagliatoNuovo].Vertices1D[1]) || (vects[Triangolo.Edges[i]].Vertices1D[1] == newopposite && vects[Triangolo.Edges[i]].Vertices1D[0] == vects[idLatoTagliatoNuovo].Vertices1D[1])) {
                Triangolo.Edges[i] = Unione.Id1D;
            }
        }


        //aggiorno matrice di adiacenza
        Matr.push_back({Triangolo.Id2D, UltimoTriangolo.Id2D});


        for (unsigned int i=0; i<3; i++) {
            if(UltimoTriangolo.Edges[i] != Unione.Id1D && UltimoTriangolo.Edges[i] != idLatoTagliatoNuovo){
                for (unsigned int j = 0; j < 2; j++) {
                    if(Matr[UltimoTriangolo.Edges[i]][j] == Triangolo.Id2D){
                        Matr[UltimoTriangolo.Edges[i]][j] = UltimoTriangolo.Id2D;
                        break;
                    }
                }

            }

        }
        for (unsigned int i=0; i<2; i++) {
            if (Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                Matr[idLatoTagliatoNuovo][i] = UltimoTriangolo.Id2D;
            }
        }


    } //fine if (stesso lato max)
    else { 


        // va fatto un if che mette l'end di latoMax = end latotagliatonuovo (caso 1)
        //                            origin latoMax = origin latotagliatovecchio (caso 2)

        for (unsigned int i = 0; i <2 ; i++) {
            if (vects[latoMax].Vertices1D[i] == vects[idLatoTagliatoVecchio].Vertices1D[0] && 1==1){
                unsigned int supporto = vects[latoMax].Vertices1D[i];
                vects[latoMax].Vertices1D[i] = vects[latoMax].Vertices1D[i-1];
                vects[latoMax].Vertices1D[i-1] = supporto;

            }
            else if (vects[latoMax].Vertices1D[i] == vects[idLatoTagliatoVecchio].Vertices1D[1] && 1==0) {
                unsigned int supporto = vects[latoMax].Vertices1D[i];
                vects[latoMax].Vertices1D[i] = vects[latoMax].Vertices1D[i+1];
                vects[latoMax].Vertices1D[i+1] = supporto;
            }
        }



        unsigned int markerMaxEdgePropa = vects[latoMax].marker1D;
        unsigned int idAltroMaxEdgePropa = 0;
        unsigned int idAltroTriPropa = 0;
        if (markerMaxEdgePropa == 0) {
            for (unsigned int i = 0; i<2; i++) {
                if (Matr[latoMax][i] != Triangolo.Id2D) {
                    idAltroTriPropa = Matr[latoMax][i];
                    idAltroMaxEdgePropa = vectt[Matr[latoMax][i]].maxedge();
                }
            }
        }



        array<unsigned int, 3> latiTriNuovoPropa = Triangolo.Edges;
        array<unsigned int, 3> vertTriNuovoPropa = Triangolo.Vertices2D;

        Vector2d midCoordPropa;
        midCoordPropa[0] = (vectp[vects[latoMax].Vertices1D[0]].Coord[0] + vectp[vects[latoMax].Vertices1D[1]].Coord[0]) *0.5;
        midCoordPropa[1] = (vectp[vects[latoMax].Vertices1D[0]].Coord[1] + vectp[vects[latoMax].Vertices1D[1]].Coord[1]) *0.5;

        unsigned int markerPPropa;
        if (vectp[vects[latoMax].Vertices1D[0]].marker0D == 0 || vectp[vects[latoMax].Vertices1D[1]].marker0D == 0) {
            markerPPropa = 0;
        }
        else {
        markerPPropa = vectp[vects[latoMax].Vertices1D[0]].marker0D; // per come sono i dati di partenza non ci sono/possono
                                                                     // essere ulteriori configurazioni
        }

        unsigned int newIndexpointPropa = vectp.size();
        Project::Cell0D newVertexPropa = Cell0D(newIndexpointPropa, markerPPropa, midCoordPropa);
        vectp.push_back(newVertexPropa);

        unsigned int oppositePropa = 0;
        for(unsigned int i = 0; i < 3; i++)
        {
            if(!(vects[latoMax].Vertices1D[0] == Triangolo.Vertices2D[i] || vects[latoMax].Vertices1D[1] == Triangolo.Vertices2D[i]))
                {
                oppositePropa = Triangolo.Vertices2D[i];
                break;
                }
        }

        vector<unsigned int> MedianaVertPropa = {oppositePropa, newVertexPropa.Id0D};

        unsigned int idNewEdgePropa = vects.size();

        unsigned int markerMedianaPropa = 0; // NON PUO' ESSERE ALTRIMENTI

        Cell1D MedianaPropa = Cell1D(idNewEdgePropa, markerMedianaPropa, MedianaVertPropa);
        vects.push_back(MedianaPropa);

        vector<unsigned int> NewSegVertPropa = {newVertexPropa.Id0D, vects[latoMax].Vertices1D[1]};


        unsigned int idNewEdgeTemp = idNewEdgePropa + 1;

        Cell1D newSegmentPropa = Cell1D(idNewEdgeTemp, vects[latoMax].marker1D, NewSegVertPropa);
        vects.push_back(newSegmentPropa);

        // aggiorno end pto medio propa
        vects[latoMax].Vertices1D[1] = newVertexPropa.Id0D; //GUARDATO FINO A QUI


        unsigned int lontano = 0;
        for (unsigned int i = 0;i<3;i++) {
            if (Triangolo.Vertices2D[i] != vects[idLatoTagliatoNuovo].Vertices1D[1] && Triangolo.Vertices2D[i] != vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                lontano = Triangolo.Vertices2D[i];
            }
        }

        // vertici effettivi nuovo propa
        for (unsigned int i = 0;i<3;i++) {
            if (vertTriNuovoPropa[i] == lontano) {
                vertTriNuovoPropa[i] = newVertexPropa.Id0D;
                break;
            }
        }

        // vertici aggiornati
        for (unsigned int i = 0;i<3;i++) {
            if (Triangolo.Vertices2D[i] != oppositePropa && Triangolo.Vertices2D[i] != lontano) {
                Triangolo.Vertices2D[i] = newVertexPropa.Id0D;
                break;
            };
        }

        // lati aggiornati
        if (vects[latoMax].Vertices1D[0] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            for (unsigned int i=0; i < 3; i++) {
                if (Triangolo.Edges[i] == latoMax){
                    Triangolo.Edges[i] = newSegmentPropa.Vertices1D[1];
                }
                if (Triangolo.Edges[i] == idLatoTagliatoVecchio){
                    Triangolo.Edges[i] = MedianaPropa.Id1D;
                }
            }

        }
        else if (vects[latoMax].Vertices1D[1] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i=0; i < 3; i++) {
                if (Triangolo.Edges[i] == idLatoTagliatoVecchio) {
                    Triangolo.Edges[i] = MedianaPropa.Id1D;
                }
            }
        }

        // lati effettivi nuovo propa
        if (vects[latoMax].Vertices1D[0] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            for (unsigned int i=0; i < 3; i++) {
                if (latiTriNuovoPropa[i] != latoMax && latiTriNuovoPropa[i] != idLatoTagliatoVecchio){
                    latiTriNuovoPropa[i] = MedianaPropa.Id1D;
                }
            }
        }
        else if (vects[latoMax].Vertices1D[1] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i=0; i < 3; i++) {
                if (latiTriNuovoPropa[i] == latoMax){
                    latiTriNuovoPropa[i] = newSegmentPropa.Id1D;
                }
                if ((vects[latiTriNuovoPropa[i]].Vertices1D[0] == oppositePropa && vects[latiTriNuovoPropa[i]].Vertices1D[1] == lontano) || (vects[latiTriNuovoPropa[i]].Vertices1D[1] == oppositePropa && vects[latiTriNuovoPropa[i]].Vertices1D[0] == lontano)) {
                    latiTriNuovoPropa[i] = MedianaPropa.Id1D;
                }
            }
        }

        unsigned int sizePenultimo = vectt.size();

        Cell2D Penultimo = Cell2D(sizePenultimo, vertTriNuovoPropa, latiTriNuovoPropa);
        vectt.push_back(Penultimo);

        //aggiorno PARZIALMENTE matrice di adiac

        if (vects[latoMax].Vertices1D[0] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
            Matr.push_back({Triangolo.Id2D, Penultimo.Id2D});
            Matr.push_back({Triangolo.Id2D, idAltroTriPropa});
            for (unsigned int i = 0; i<2; i++) {
                if (Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                    Matr[idLatoTagliatoNuovo][i] = Penultimo.Id2D;
                }
            }
        }
        else if (vects[latoMax].Vertices1D[1] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            Matr.push_back({Triangolo.Id2D, Penultimo.Id2D});
            for (unsigned int i = 0; i<2; i++) {
                if (Matr[idLatoTagliatoVecchio][i] == Triangolo.Id2D) {
                    Matr[idLatoTagliatoVecchio][i] = Penultimo.Id2D;
                }

            }
        }


        // collego i due punti medi

        array<unsigned int, 3> latiTriResiduoPropa = Penultimo.Edges;
        array<unsigned int, 3> vertTriResiduoPropa = Penultimo.Vertices2D;


        unsigned int sizeUnione = vects.size();
        unsigned int markerUnione = 0;
        vector<unsigned int> vettoreUnione = {vects[idLatoTagliatoVecchio].Vertices1D[1], newVertexPropa.Id0D};

        Cell1D Unione = Cell1D(sizeUnione, markerUnione, vettoreUnione);
        vects.push_back(Unione);

        if (vects[latoMax].Vertices1D[0] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {

            for (unsigned int i = 0; i<3; i++) {
                // aggiorno vertici
                if (Penultimo.Vertices2D[i] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                    Penultimo.Vertices2D[i] = vects[idLatoTagliatoVecchio].Vertices1D[1];
                }
                if (vertTriResiduoPropa[i] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
                    vertTriResiduoPropa[i] = vects[idLatoTagliatoNuovo].Vertices1D[0];
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

            unsigned int sizeResiduoT = vectt.size();

            Cell2D Residuo = Cell2D(sizeResiduoT, vertTriResiduoPropa, latiTriResiduoPropa);
            vectt.push_back(Residuo);

            // aggiorno matrice adiacenza
            Matr.push_back({Penultimo.Id2D, Residuo.Id2D});
            for (unsigned int i = 0; i<2; i++) {
                if (Matr[latoMax][i] == Triangolo.Id2D) {
                    Matr[latoMax][i] = Residuo.Id2D;
                }
                if (Matr[idLatoTagliatoVecchio][i] == Triangolo.Id2D) {
                    Matr[idLatoTagliatoVecchio][i] = Residuo.Id2D;
                }
            }
        }
        else if (vects[latoMax].Vertices1D[1] == vects[idLatoTagliatoNuovo].Vertices1D[1]) {
            for (unsigned int i = 0; i<3; i++) {
                if (Penultimo.Vertices2D[i] == newSegmentPropa.Vertices1D[1]) {
                    Penultimo.Vertices2D[i] = vects[idLatoTagliatoVecchio].Vertices1D[1];
                }
                if (vertTriResiduoPropa[i] == vects[idLatoTagliatoVecchio].Vertices1D[0]) {
                    vertTriResiduoPropa[i] = vects[idLatoTagliatoNuovo].Vertices1D[0];
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

            unsigned int sizeResiduo = vectt.size();

            Cell2D Residuo = Cell2D(sizeResiduo, vertTriResiduoPropa, latiTriResiduoPropa);
            vectt.push_back(Residuo);

            // aggiorno matrice adiacenza
            Matr.push_back({Residuo.Id2D, idAltroTriPropa});
            Matr.push_back({Penultimo.Id2D, Residuo.Id2D});
            for (unsigned int i = 0; i<2; i++) {
                if (Matr[idLatoTagliatoNuovo][i] == Triangolo.Id2D) {
                    Matr[idLatoTagliatoNuovo][i] = Residuo.Id2D;
                }
            }
        }

        // void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Cell2D Triangolo, unsigned int latoMax){
        // in questo caso è ricorsiva
        if (markerMaxEdgePropa == 0) {
            Cell2D AltroTriPropa = vectt[idAltroTriPropa];
            Propagazione(latoMax, newSegmentPropa.Id1D, AltroTriPropa, idAltroMaxEdgePropa, vectp, vects, vectt, Matr);
        }
    } // fine else (lato lungo diverso dal precedente)
} // fine Propagazione
//Export
bool TriangularMesh::ExportCell0Ds(TriangularMesh& mesh){
        ofstream file;
        string outputFile= "ExportCell0Ds.csv";
        file.open(outputFile);
        if (!outputFile.is_open()) {
            cerr << "Failed to open the file." << endl;
            return false;
        }
        file << "id marker coord" <<endl; 

        for(unsigned int i=0; i<mesh.numbercell0D; i++)
        {
            outputFile << mesh.Cell0DId[i] << " " << mesh.Cell0DMarkers.find(marker) << " " ;
            outputFile << mesh.Cell0DCoordinates[mesh.Cell0DId[i]][0] << " " << mesh.Cell0DCoordinates[mesh.Cell0DId[i]][1] << "\n";
        }

        file.close();

        return true;
    }
 bool TriangularMesh::ExportCell1Ds(TriangularMesh& mesh){
        ofstream file;
        string outputFile = "ExportCell1Ds.csv";
        file.open(outputFile);
        if (!outputFile.is_open()) {
            cerr << "Failed to open the file." << endl;
            return false;
        }
        file << "Id Marker Vertices" <<endl; 

        for(unsigned int i=0; i<mesh.numbercell1D; i++)
        {
            outputFile << mesh.Cell1DId[i] << " " << mesh.Cell1DMarkers.find(marker) << " ";
            outputFile << mesh.Cell1DVertices[mesh.Cell1DId[i]][0] << " " << mesh.Cell1DVertices[mesh.Cell1DId[i]][1] << "\n";
        }

        file.close();

        return true;
    }
bool TriangularMesh::ExportCell2Ds(TriangularMesh& mesh){
    ofstream file;
    string outputFile= "ExportCell2Ds.csv";
    file.open(outputFile);
        if (!outputFile.is_open()) {
            cerr << "Failed to open the file." << endl;
            return false;
        }
        file << "Id Vertices Edges" <<endl; 

        for(unsigned int i=0; i<mesh.numbercell2D; i++)
        {
            outputFile << mesh.Cell2DId[i] << " ";

            outputFile << mesh.Cell2DVertices2D[mesh.Cell2DId[i]][0] << " " << mesh.Cell2DVertices2D[mesh.Cell2DId[i]][1] << " " << mesh.Cell2DVertices2D[mesh.Cell2DId[i]][2] << " ";
            outputFile << mesh.Cell2DEdges[mesh.Cell2DId[i]][0] << " " << mesh.Cell2DEdges[mesh.Cell2DId[i]][1] << " " << mesh.Cell2DEdges2D[mesh.Cell2DId[i]][2] <<endl;
        }

        file.close();

        return true;
}
//}


 // fine namespace Cells


//-----------------------------------------------------------------------
