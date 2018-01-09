#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#include "matlab/mat.h"
#include "matlab/mex.h"
#include "matlab/matrix.h"

#pragma comment(lib, "lib/libmex.lib")
#pragma comment(lib, "lib/libmx.lib")
#pragma comment(lib, "lib/libmat.lib")

struct point3d
{
  float x, y, z;
  point3d() {}
  point3d(float p[3]) { x = p[0], y = p[1], z = p[2]; }
};

/// A node element, which contains who is its father, and the list of point to reach it
class CTreeNode
{
public:
  std::vector<point3d> m_vecPoints;
  int m_iParent;
  CTreeNode() : m_iParent { -1 } {  }
};


class CArrayTree
{
private:
  std::vector<CTreeNode> m_tree;
  /// Determines the equivalent subscript values corresponding to the absolute index dimension of a multidimensional array.
  /// @param siz size of the N-dimensional matrix
  /// @param N the dimensions of the matrix
  /// @param idx index in linear format
  /// @param sub the output - subscript values written into an N - dimensional integer array (created by caller)
  void ind2sub(const size_t *siz, int N, int idx, int *sub)
  {
    // Calculate the cumulative product of the size vector
    std::vector<int> cumProd(N);
    int nextProd = 1;
    for (int k = 0; k < N; k++)
    {
      cumProd[k] = nextProd;
      nextProd *= siz[k];
    }

    // Validate range of linear index
    int maxInd = cumProd[N - 1] * siz[N - 1] - 1;	//assuming that fit into a int data type
    if ((idx < 0) || (idx > maxInd))
    {
      std::cout << "bad linear index maxsum::ind2sub" << std::endl;	//some kind of message
      return;
    }

    //the computation
    int remainder = idx;
    for (int k = N - 1; k >= 0; k--)
    {
      sub[k] = remainder / cumProd[k];
      remainder = remainder % cumProd[k];
    }
  }
public:
  CArrayTree() {  }
  ~CArrayTree() {  }

  // Load a .mat file into the tree structure
  // @param name of the filename to open
  bool loadMATFile(string strFileName)
  {
    MATFile* m_pFileMat = matOpen(strFileName.c_str(), "r");
    if (m_pFileMat == nullptr) return false;

    //load variable G
    mxArray* m_pArray = matGetVariable(m_pFileMat, "GFinalDir");
    //get size volume
    //const mwSize *dimsE = mxGetDimensions(pArrToE);
    mxArray *pArraySkeleton = matGetVariable(m_pFileMat, "skel");
    const mwSize *dimskel = mxGetDimensions(pArraySkeleton);

    //get the field e
    mxArray *pArrToE = mxGetField(m_pArray, 0, "e");
    const mwSize *dimsE = mxGetDimensions(pArrToE);
    int iWidth = dimsE[0];
    int iHeight = dimsE[1];
    double* pValue = mxGetPr(pArrToE);

    //get the field we
    mxArray *pArrToWE = mxGetField(m_pArray, 0, "we");  //this is a cell

    //create the tree
    m_tree.resize(iWidth + 1);

    //create matrix E entirely 
    //matlab stores the data column-order - IMPORTANT! -
    std::vector<std::vector<int>> m_matrixE;
    m_matrixE.resize(iHeight);
    for (mwIndex y = 0; y < iHeight; y++)
      m_matrixE[y].resize(iWidth);
    
    //start the loop
    for (mwIndex x = 0; x < iWidth; x++)
    {
      for (mwIndex y = 0; y < iHeight; y++)
      {
        m_matrixE[x][y] = pValue[(y*iWidth) + x];
        if (m_matrixE[x][y] >= 1)
        {
          m_tree[(y+1)].m_iParent = (x + 1);  //copy the value

          //now, extract points
          mxArray *pCellContent = mxGetCell(pArrToWE, [](int x, int y, int width) -> size_t { //return the linear index of x, y
            return (width*y) + x;
          }(x, y, iWidth));
          //now, we have the number of elements of GFinalDir.we
          const mwSize pNumElements = mxGetDimensions(pCellContent)[1];  //always is 1 x k, where k is the size. then, the position 0 is always 1
          double* pValueWE = mxGetPr(pCellContent);
          for (int k = 0; k < pNumElements; k++) 
          {
            int point[3];
            point3d thepoint;
            ind2sub(dimskel, 3, pValueWE[k], point);
            thepoint.x = point[0];
            thepoint.y = point[1];
            thepoint.z = point[2];
            m_tree[(y+1)].m_vecPoints.push_back(thepoint);
          }
        }
      }
    }
    return true;
  }

  // From root (node 1), this prints all points to reach until iNode
  // @param iNode destination node
  void printPointsPath(int iNode) 
  {
    int iTempNode = iNode;
    std::vector<int> path;
    while (iTempNode != 1) //until reach root
    {
      path.push_back(iTempNode);
      iTempNode = m_tree[iTempNode].m_iParent;
    }
    path.push_back(1);
    std::reverse(path.begin(), path.end());

    for_each(path.begin(), path.end(), [=](int index) 
    //for_each(path.begin(), path.begin() + 5, [=](int index)
    {
      for_each(m_tree[index].m_vecPoints.begin(), m_tree[index].m_vecPoints.end(), [](point3d p) {
        cout << p.x << " " << p.y << " " << p.z << endl;
      });
    });
  }

  // From root (node 1), prints the sequence in format a --> b --> ... --> iNode
  void printPath(int iNode)
  {
    int iTempNode = iNode;
    std::vector<int> path;
    while(iTempNode != 1) //until reach root
    {
      path.push_back(iTempNode);
      iTempNode = m_tree[iTempNode].m_iParent;
    }
    path.push_back(1);
    std::reverse(path.begin(), path.end()); //are stored in backward order, then reverse it

    for_each(path.begin(), path.end(), [](int index) {
      cout << index << " --> ";
    });
  }
};

int main() 
{
  string strFilename = "C:\\Users\\eramirez\\Desktop\\CARS\\Cases\\28_INSP_CPAP_Graph_skel_vol.mat";
  CArrayTree app;
  if (app.loadMATFile(strFilename)) 
  {
    //app.printPath(2);
    app.printPointsPath(40);
  }
  else 
  {
    cout << "some error reading" << endl;
  }
  return EXIT_SUCCESS;
}