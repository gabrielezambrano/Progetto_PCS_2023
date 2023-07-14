#ifndef __EMPTY_H
#define __EMPTY_H

#include <iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <algorithm>
#include "cmath"

using namespace std;
using namespace Eigen;

namespace Project {




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
            Cell1D(unsigned int id, unsigned int marker, vector<unsigned int> vertices);
            double LengthEdge();
 };


    class Cell2D{
        public:
            unsigned int LengthEdges;
            unsigned int Id2D;
            array<unsigned int, 3> Vertices2D;
            array<unsigned int, 3> Edges;
            Cell2D(unsigned int id, array<unsigned int, 3> Vertices2D, array<unsigned int, 3> Edges);
            unsigned int maxedge();
            double Area();
        };

  class TriangularMesh{
    public:
        unsigned int numbercell0D;
        vector<Project::Cell0D> vectp = {};



        unsigned int numbercell1D;
        vector<Project::Cell1D> vects = {};
        vector<double> LengthEdges = {};


        unsigned int numbercell2D;
        std::vector<vector<unsigned int>> LenghtMax = {};
        vector<Project::Cell2D> vectt = {};

        TriangularMesh() = default;
    };

    class MatrAdiac{
    public:
        vector<vector<unsigned int>> Matr;
        MatrAdiac(vector<Project::Cell2D> vectt, vector<Project::Cell1D> vects);
    };

 }
   
    bool ImportCell0Ds();
    bool ImportCell1Ds();
    bool ImportCell2Ds();
    bool ExportCell0Ds();
    bool ExportCell1Ds();

  void Bisect(Project::Cell2D triangleToBisect, vector<Project::Cell0D> vectp, vector<Project::Cell1D> vects, vector<Project::Cell2D> vectt, vector<vector<unsigned int>> Matr);
  void Propagazione(unsigned int idLatoTagliatoVecchio, unsigned int idLatoTagliatoNuovo, Project::Cell2D Triangolo, unsigned int latoMax, vector<Project::Cell0D> vectp, vector<Project::Cell1D> vects, vector<Project::Cell2D> vectt, vector<vector<unsigned int>> Matr);

    // tolleranze
    constexpr double max_tol(const double& x, const double& y)
        {
        return x > y ? x : y;
        }

    static constexpr double tol1D = 1.0e-12;
    static constexpr double tol2D = max_tol(tol1D * tol1D, numeric_limits<double>::epsilon());


    inline bool operator<(Project::Cell2D& t1,  Project::Cell2D& t2)
    {
      return t1.Area() < t2.Area() + tol2D * max(t1.Area(), t2.Area());
    }

    inline bool operator>=(Project::Cell2D& t1,Project::Cell2D& t2)
    {
      return !(t1 < t2);
    }



 


#endif // __EMPTY_H
