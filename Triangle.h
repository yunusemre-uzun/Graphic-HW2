#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <vector>

class Triangle
{
public:
    int vertexIds[3];
    bool isCulled = false;
    std::vector<int> clipped_starting_point_indexes;
    std::vector<int> clipped_ending_point_indexes;

    Triangle();
    Triangle(int vid1, int vid2, int vid3);
    Triangle(const Triangle &other);

    int getFirstVertexId();
    int getSecondVertexId();
    int getThirdVertexId();

    void setFirstVertexId(int vid);
    void setSecondVertexId(int vid);
    void setThirdVertexId(int vid);
};


#endif