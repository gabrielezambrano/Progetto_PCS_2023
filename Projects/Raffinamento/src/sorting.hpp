ifndef SORTING_HPP
#define SORTING_HPP


#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;

namespace Sorting
{

template <typename T>
void MakeHeap(vector<T>& vecttSupp, int i){



    int max = i;
    unsigned int l = 2 * i + 1;
    unsigned int r = 2 * i + 2;

    if (l < vecttSupp.size() && vecttSupp[l] < vecttSupp[max])
    {
        max = l;
    }

    if (r < vecttSupp.size() && vecttSupp[r] < vecttSupp[max])
    {
        max = r;
    }

    if (max != i)
    {
        swap(vecttSupp[i], vecttSupp[l]);
        MakeHeap(vecttSupp, i);
    }
}

template<typename T>
void HeapSort(vector<T>& vecttSupp, vector<T>& vectt){

    vecttSupp.clear();
    for (unsigned int k = 0; k < vectt.size(); k++){
        vecttSupp.push_back(vectt[k]);
    }

    for (int i = vecttSupp.size() / 2 - 1; i >= 0; i--)
    {
        MakeHeap(vecttSupp, i);
    }
    for (int i = vecttSupp.size() - 1; i >= 0; i--)
    {
        swap(vecttSupp[0], vecttSupp[i]);
        MakeHeap(vecttSupp, 0);
    }
}
}



#endif // SORTING_HPP


