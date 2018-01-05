#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "matlab/mat.h"
#include "matlab/mex.h"
#include "matlab/matrix.h"

#pragma comment(lib, "lib/libmex.lib")
#pragma comment(lib, "lib/libmx.lib")
#pragma comment(lib, "lib/libmat.lib")

class Mat2CPP 
{
private:
  std::vector<std::vector<int>> m_matrixE;
  MATFile* m_pFileMat; 
  mxArray *m_pArray;
public:
  Mat2CPP()  
  {
    m_pFileMat = nullptr;
    m_pArray = nullptr;
  }
  
  bool loadMATFile(string strFileName)
  {
    m_pFileMat = matOpen(strFileName.c_str(), "r");
    if (m_pFileMat == nullptr) return false;

    //load variable G
    m_pArray = matGetVariable(m_pFileMat, "G");
    //get the field e
    mxArray *pArrToE = mxGetField(m_pArray, 0, "e");
    const mwSize *dimsE = mxGetDimensions(pArrToE);
    int iWidth = dimsE[0];
    int iHeight = dimsE[1];
    double* pValue = mxGetPr(pArrToE);
    //create matrix E
    m_matrixE.resize(iHeight);
    for (mwIndex y = 0; y < iHeight; y++)
    {
      m_matrixE[y].resize(iWidth);
      for (mwIndex x = 0; x < iWidth; x++) 
      {
        if (pValue[(y*iWidth) + x] > 0) 
        {
          cout << "(" << x << "," << y << ") = " << pValue[(y*iWidth) + x] << endl;
        }
      }
    }
    return true;
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
  string strFilename = "C:\\Users\\esmitt\\Desktop\\CARS\\Cases\\14_INSP_CPAP_Graph_skel_vol.mat";
  Mat2CPP app;
  if (app.loadMATFile(strFilename)) 
  {
  
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