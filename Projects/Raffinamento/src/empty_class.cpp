#include "empty_class.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include <fstream>

using namespace std;

vector<Cells::Cell0D> vectp = {};
vector<Cells::Cell1D> vects = {};
vector<Cells::Cell2D> vectt = {};
Cells::TriangularMesh mesh;

bool ImportCell0Ds();
bool ImportCell1Ds();
bool ImportCell2Ds();

namespace Cells
{



//Cell0D::Cell0D(const unsigned int& id,
//               const Vector2d& coordinates0D)
Cell0D::Cell0D(unsigned int id, unsigned int marker, Vector2d coord)
    {
    unsigned int marker0D = marker;
    unsigned int Id0D = id;
    Vector2d Coord = coord;
    };



Cell1D::Cell1D(unsigned int id,unsigned int marker,Vector2i vertices)
    {
    unsigned int Id1d = id;
    unsigned int marker1D = marker;
    Vector2i Vertices1d = vertices;
    };

    
//metterei double anziché void
double Cells::Cell1D::LengthEdge(){
    Vector2d coordOrigin = mesh.coordinates0D[Vertices1D[0]];
    Vector2d coordEnd= mesh.coordinates0D[Vertices1D[1]];
    //LengthEdges = (coordEnd-coordOrigin).norm();
    double len = sqrt(pow(coordOrigin[0]-coordEnd[0], 2)+pow(coordOrigin[1] - coordEnd[1], 2));
    return len;
    }




//metterei double anziché void
//double Cell1D::LengthEdge(){
//    Vector2d coordOrigin = {0, 0};
//    Vector2d coordEnd= {0, 0};
//    bool a = false;
//    bool b = false;
//    for(unsigned int i = 0;i<mesh.numbercell0D;i++)
//       {
//        if(Vertices1D[0]==mesh.id0D[i])
//            {
//            coordOrigin = mesh.coordinates0D[i];
//            a = true;
//          break;
//            }
//        if(Vertices1D[1]==mesh.id0D[i])
//            {
//            coordEnd = mesh.coordinates0D[i];
//            b = true;
//          break;
//            }
//        if(a == true && b == true)
//        {
//            break;
//        }
//        //LengthEdges = (coordEnd-coordOrigin).norm();
//        double len = sqrt(pow(coordOrigin[0]-coordEnd[0], 2)+pow(coordOrigin[1] - coordEnd[1], 2));
//        return len;
//        }
//}    
    
  Cell2D::Cell2D(unsigned int id,array<unsigned int, 3> Vertices, array<unsigned int, 3> Edges)
    {
    unsigned int Id2D = id;
    array<unsigned int, 3> Vertices2D = Vertices;
    array<unsigned int, 3> Edges2D = Edges;
    };
 
    
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

//unsigned int Cell2D::maxedge(){
//    unsigned int indmax = 0;
//    double max = 0.0;
//    array<double, 3> lenedges = {0, 0, 0};
//    bool a = false;
//    bool b = false;
//    bool c = false;
//    for(unsigned int i = 0;i<mesh.numbercell1D;i++)
//       {
//        if(Edges[0]==mesh.id1D[i])
//            {
//            lenedges[0] = mesh.LengthEdges[i];
//            a = true;
//          break;
//            }
//        if(Edges[1]==mesh.id1D[i])
//            {
//            lenedges[1] = mesh.LengthEdges[i];
//            b = true;
//          break;
//            }
//        if(Edges[2]==mesh.id1D[i])
//            {
//            lenedges[2] = mesh.LengthEdges[i];
//            c = true;
//          break;
//            }
//        if(a == true && b == true && c == true)
//        {
//            break;
//        }
//    }
//    for (unsigned int i = 0; i<3; i++)
//    {
//        if(lenedges[i] > max)
//        {
//            max = lenedges[i];
//            indmax = Edges[i]; //NON SO SE HA PIU' SENSO CHE RESTITUISCA L'ID DEL LATO O IL SUO INDICE NEL VETTORE DA TRE
//        }
//    }
//    return indmax;

//}    
    
  double Cell2D::Area(){
             //Formula dell'area di Gauss
             A_12 = (x_1*y_2) - (y_1*x_2);
             A_23 = (x_2*y_3) - (y_2*x_3);
             A_31 = (x_3*y_1) - (y_3*x_1);
             Area = abs((A_12+A_23+A_31)/2);
  }
    
    
//vector<vector<unsigned int>> MatrAdiac(mesh.numbercell1D);
//bool Raffinamento::CreationMatrAdiac() {
//    for (unsigned int i = 0; i < mesh.numbercell2D; i++) {
//        for (unsigned int j = 0; j < 3; j++) {
//            MatrAdiac[mesh.vectt[i].Edges[j]].push_back(mesh.vectt[i].Id2D);
//        }
//    }
//return true;
//}    
    
//-----------------------------------------------------------------
    
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

  unsigned int NumberCell0D = listLines.size();
  unsigned int numberCell0D = NumberCell0D;

  if (NumberCell0D == 0)
  {
    cerr << "There is no cell 0D" << endl;
    return false;
  }

  mesh.id0D.reserve(numberCell0D);
  mesh.coordinates0D.reserve(numberCell0D);

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
    vectp.push_back(point);

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

  unsigned int NumberCell1D  = listLines.size();
  unsigned int numberCell1D = NumberCell1D;

  if (NumberCell1D == 0)
  {
    cerr << "There is no cell 1D" << endl;
    return false;
  }

  mesh.id1D.reserve(numberCell1D);
  mesh.vertices1D.reserve(numberCell1D);

  for (const string& line : listLines)
  {
    istringstream converter(line);

    unsigned int id;
    unsigned int marker;
    Vector2i vertices;

    converter >>  id >> marker >> vertices(0) >> vertices(1);
    Cells::Cell1D segment = Cells::Cell1D(id,marker,vertices);
    vects.push_back(segment);
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

  unsigned int NumberCell2D = listLines.size();
  unsigned int numberCell2D = NumberCell2D;

  if (NumberCell2D == 0)
  {
    cerr << "There is no cell 2D" << endl;
    return false;
  }

  mesh.id2D.reserve(NumberCell2D);
  mesh.vertices2D.reserve(NumberCell2D);
  mesh.edges2D.reserve(NumberCell2D);

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
    vectt.push_back(triangle);
    mesh.id2D.push_back(id);
    mesh.vertices2D.push_back(vertices);
    mesh.edges2D.push_back(edges);
  }
  file.close();
  return true;
}

