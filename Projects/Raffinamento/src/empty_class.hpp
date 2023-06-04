//{
// class Empty
//  {
//    public:
//      void Show() const { std::cout<< "Hello world;"<< std::endl; }
//  };
//}
//
//#endif // __EMPTY_H


//#ifndef __EMPTY_H
//#define __EMPTY_H

#include <iostream>
#include "Eigen/Eigen"
#include <fstream>

using namespace std;
using namespace Eigen;

namespace Cells {

class TriangularMesh{
public:
    unsigned int numbercell0D;
    vector<unsigned int> id0D;
    vector<Vector2d> coordinates0D;


    unsigned int numbercell1D;
    vector<unsigned int> id1D;
    vector<Vector2i> vertices1D;
    Eigen::VectorXd LengthEdges = {};


    unsigned int numbercell2D;
    vector<unsigned int> id2D;
    vector<array<unsigned int, 3>> edges2D;
    vector<array<unsigned int, 3>> vertices2D;
    std::vector<vector<unsigned int>> LenghtMax = {};


//class Raffinamento{
//    public:
//        vector<vector<unsigned int>> MatrAdiac;
//        bool CreationMatrAdiac();
//    };    
    
};
 class Cell0D {
    private:
        unsigned int marker0D;

    public:
        unsigned int Id0D;
        Vector2d Coord;
    Cell0D(unsigned int id,unsigned int marker,Vector2d coord);
    };


    class Cell1D {
    private:
        unsigned int marker1D;

    public:
        unsigned int Id1D;
        vector<unsigned int> Vertices1D;
        Cell1D(unsigned int id, unsigned int marker, Vector2i vertices);
        double LengthEdge();
 };


    class Cell2D{
    public:
        unsigned int LengthEdges;
        unsigned int Id2D;
        array<unsigned int, 3> Vertices2D;
        array<unsigned int, 3> Edges;
        Cell2D(unsigned int id,array<unsigned int, 3> Vertices, array<unsigned int, 3> Edges);
        double MaxEdge();
        double Area();
    };

 }

