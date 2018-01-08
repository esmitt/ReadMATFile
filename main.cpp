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

//class CDataTree
//{
//public:
//  int m_iMatLabIndex;
//  int m_inPoints;
//  CDataTree() 
//  {
//    m_iMatLabIndex = -1;
//    m_inPoints = 0;
//  }
//};
//
//class CMatTreeNode
//{
//public:
//  CDataTree m_data;
//  int m_iIndex;
//  int m_iParent;
//  int m_iChildL;
//  int m_iChildR;
//
//  CMatTreeNode()
//  {
//    m_iIndex = -1;
//    m_iParent = -1;
//    m_iChildL = -1;
//    m_iChildR = -1;
//  }
//
//  void printNode()
//  {
//    cout << "index: " << m_iIndex << endl;
//    cout << "Matlab-index: " << m_data.m_iMatLabIndex << endl;
//    cout << "points: " << m_data.m_inPoints << endl;
//    cout << "parent " << m_iParent << endl;
//    cout << "LChild-index " << m_iChildL << endl;
//    cout << "RChild-index " << m_iChildR << endl;
//  }
//};
//
//enum class NODE_TYPE {ROOT, LEFT, RIGHT, NONE};
typedef struct point3d
{
  float x, y, z;
  point3d() {}
  point3d(float p[3]) { x = p[0], y = p[1], z = p[2]; }
};

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

class CTreeNode
{
public:
  std::vector<point3d> m_vecPoints;
  int m_iParent;
  CTreeNode() : m_iParent { -1 } {  }
};

class Mat2CPP 
{
private:
  std::vector<std::vector<int>> m_matrixE;
  MATFile* m_pFileMat;
  mxArray *m_pArray;

//  std::vector<CMatTreeNode> vecTree;
  std::vector<CTreeNode> m_tree;
public:
  Mat2CPP()  
  {
    m_pFileMat = nullptr;
    m_pArray = nullptr;
  }
  
  //void printTree() 
  //{
  //  for_each(vecTree.begin(), vecTree.end(), [](CMatTreeNode& node) {
  //    if(node.m_iIndex > 0)
  //      node.printNode();
  //  });
  //}

  bool loadMATFile(string strFileName)
  {
    m_pFileMat = matOpen(strFileName.c_str(), "r");
    if (m_pFileMat == nullptr) return false;

    //load variable G
    m_pArray = matGetVariable(m_pFileMat, "GFinalDir");
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
    mxArray *pArrToWE = mxGetField(m_pArray, 0, "we");

    //create the tree
    //vecTree.resize(iWidth);
    m_tree.resize(iWidth + 1);

    //create matrix E entirely 
    //matlab stores the data column-order
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
          //m_tree[(y+1)].m_inPoints = m_matrixE[x][y];
          m_tree[(y+1)].m_iParent = (x + 1);
          //now, extract points
          mxArray *pCellContent = mxGetCell(pArrToWE, [](int x, int y, int width) -> size_t { //return the linear index of x, y
            return (width*y) + x;
          }(x, y, iWidth));
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
    //printTree();
    return true;
  }

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

    //for_each(path.begin(), path.end(), [=](int index) 
    for_each(path.begin(), path.begin() + 2, [=](int index)
    {
      for_each(m_tree[index].m_vecPoints.begin(), m_tree[index].m_vecPoints.end(), [](point3d p) {
        cout << p.x << " " << p.y << " " << p.z << endl;
      });
    });
  }

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
    std::reverse(path.begin(), path.end());
    for_each(path.begin(), path.end(), [](int index) {
      cout << index << " --> ";
    });
  }

  void printPath(int iStartNode, int iEndNode)
  {
    if (iEndNode != iStartNode) 
    {
      printPath(m_tree[iStartNode].m_iParent, iEndNode);
      cout << iStartNode << " --> ";
    }
  }

  //void compactTree() 
  //{
  //  for (int iIndex = 1; iIndex < vecTree.size(); iIndex++) 
  //  {
  //    if(vecTree[iIndex].m_iIndex > -1) //if exists
  //    {
  //      int iRoot = vecTree[iIndex].m_iIndex;
  //      //check if the root is child on another subtree
  //      for (int child = 1; child < vecTree.size(); child++) 
  //      {
  //        if (child != iIndex && vecTree[child].m_iIndex > -1) //if exists
  //        {
  //          //check if exist on left child
  //          if (iRoot == vecTree[child].m_iChildL) 
  //          {
  //            vecTree[iIndex].m_iParent = child;
  //            break;
  //          }
  //          if (iRoot == vecTree[child].m_iChildR)
  //          {
  //            vecTree[iIndex].m_iParent = child;
  //            break;
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  void print() 
  {
    for (mwIndex x = 0; x < m_matrixE.size(); x++)
      for (mwIndex y = 0; y < m_matrixE.size(); y++)
        if(m_matrixE[x][y] >= 1)
          cout << "(" << (x + 1) << "," << (y + 1) << ") = " << m_matrixE[x][y] << endl;
  }

  ~Mat2CPP() 
  {
    if(m_pFileMat != nullptr)
      matClose(m_pFileMat);
    if(m_pArray != nullptr)
      mxDestroyArray(m_pArray);
  }
};

int main() 
{
  string strFilename = "C:\\Users\\eramirez\\Desktop\\CARS\\Cases\\28_INSP_CPAP_Graph_skel_vol.mat";
  Mat2CPP app;
  if (app.loadMATFile(strFilename)) 
  {
    //app.print();
    //app.printPath(2);
    app.printPointsPath(2);
  }
  else 
  {
    cout << "some error reading" << endl;
  }
  //MATFile* m_pFileMat = matOpen(strFilename.c_str(), "r");
  //if (m_pFileMat == nullptr)
  //  return EXIT_FAILURE;
  //
  //mxArray *pArr = matGetVariable(m_pFileMat, "G");
  //mwSize num = mxGetNumberOfFields(pArr);	//number of fields

  //mxArray *pArrToE = mxGetField(pArr, 0, "e");	//get the field 
  ////mwSize iDimE = mxGetNumberOfElements(pArrToE);
  //const mwSize *dims = mxGetDimensions(pArrToE);
  //double* pValue = mxGetPr(pArrToE);
  //for (mwIndex j = 0; j < dims[1]; j++)
  //  for (mwIndex i = 0; i < dims[0]; i++) 
  //  {
  //      //cout << (i + 1) << " " << (j + 1) << endl;
  //    if (pValue[(j*dims[0]) + i] > 0) 
  //    {
  //      cout << "(" << i << "," << j << ") = " << pValue[(j*dims[0]) + i] << endl;
  //    }
  //  }
  //double *pr = mxGetPr(pArrToE);

  //mxArray *pArrToWE = mxGetField(pArr, 0, "we");
  //mxArray *cellElement;
  //const mwSize *dims = mxGetDimensions(pArrToWE);
  //
  //for (mwIndex jcell = 0; jcell < dims[0]; jcell++)
  //{
  //  cellElement = mxGetCell(pArrToWE, jcell);
  //  if (cellElement != nullptr)
  //  {
  //    double* p = mxGetPr(cellElement);
  //    mexPrintf("The content at is %g\n", *p);
  //  }
  //}

//  mxDestroyArray(pArr);
//  matClose(m_pFileMat);
  return EXIT_SUCCESS;
}