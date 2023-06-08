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


    class Cell0D {


        public:
            unsigned int marker0D;
            unsigned int Id0D;
            Vector2d Coord;
        Cell0D(unsigned int id, unsigned int marker, Vector2d coord);
        };


    class Cell1D {

        public:
            unsigned int marker1D;
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
            unsigned int maxedge();
            double Area();
        };

    class TriangularMesh{
    public:
        unsigned int numbercell0D;
        vector<Cells::Cell0D> vectp = {};



        unsigned int numbercell1D;
        vector<Cells::Cell1D> vects = {};
        vector<double> LengthEdges = {};


        unsigned int numbercell2D;
        std::vector<vector<unsigned int>> LenghtMax = {};
        vector<Cells::Cell2D> vectt = {};

    };

    class MatrAdiac{
    public:
        vector<vector<unsigned int>> Matr;
        MatrAdiac();
    };

 }
