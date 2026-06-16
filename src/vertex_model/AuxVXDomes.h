/*
    *******************************************************************************
    Copyright (c) 2017-2025
    Authors/Contributors: Adam Ouzeri, Sohan Kale, Alejandro Torres-Sánchez, Daniel Santos-Oliván
    *******************************************************************************
    This file is part of agvm_domes
    Project homepage: https://github.com/aouzeri
    Distributed under the GNU General Public License, see the accompanying
    file LICENSE or https://opensource.org/license/gpl-3-0.
    *******************************************************************************
*/

#ifndef AUXVERTEXMODEL_H
#define AUXVERTEXMODEL_H

#include <iostream>
#include <mpi.h>


#include <unistd.h> //FIXME sleep function

//#include "DistributedClass.h"
//#include "TypeDefs.h"

#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_HiPerProblem.h"
#include "hl_DistributedClass.h"
#include "hl_HiperIndexConverter.h"

#include "hl_Math.h"
#include "hl_Array.h"
#include "hl_Geometry.h"
#include "hl_Timer.h"

#include <string>


class faceMesh
{
    public:
        int faceType = 0;
        int nDim     = 3;
        int nVert    = 3;
        int nPts     = {};
        int nElem    = {};
        int DOF      = 3;
        int gDOF     = 6; // G values, pressure, density, luminal pressure
        int gAux     = 6;

        std::vector<int>    nBorFaces     = {};
        std::vector<int>    connec        = {};
        std::vector<double> x_nodes       = {};
        std::vector<int>    idN_creases   = {};
        std::vector<std::vector< std::pair<int,int>>> masterNodes = {};
        std::vector<std::vector< double >> masterNodesPerFlag = {};
        std::vector<bool> cn_nodes = {};
        double centroidPositionInFaceMesh [3] = {};


        Teuchos::RCP<hiperlife::FieldStruct> nodeDOFs = Teuchos::null;
        Teuchos::RCP<hiperlife::FieldStruct> auxF     = Teuchos::null;
        Teuchos::RCP<hiperlife::FieldStruct> gDOFS    = Teuchos::null;
        Teuchos::RCP<hiperlife::FieldStruct> auxV     = Teuchos::null;

        std::vector<double> u_nodes;
        std::vector<double> globalU;
        std::vector<double> globalV;

        void loadMesh( std::string fileName);

};

class tissueMesh: public hiperlife::DistributedData
{

    public:

        int setConstraintflag        = 1;
        int nCells                   = 0;
        int nFaces                   = 0;
        std::vector<faceMesh> meshes = {};
        std::vector< int> nPointsinFace = {};

        // FIXME: parallelise
        std::vector<vector<int >> facesInCell = {};
        std::vector<vector<int >> NeighboursCells = {};

		void loadMesh( std::string fileName);
        void readNeighbouringCells(const string fileName);
        void readfaceIDsinCells(const string fileName);
        void calculateDistanceBetweenCellCentroids(Teuchos::RCP<hiperlife::ParamStructure> paramStr,  Teuchos::RCP<hiperlife::HiPerProblem> hiperProbl);
        void setNumberCells(const int numberCells);
	    void setNPointsInCell(Teuchos::RCP<hiperlife::ParamStructure> paramStr);
        void updateFacesCentroids(Teuchos::RCP<hiperlife::ParamStructure> paramStr);

    faceMesh* operator [] (int i)
        {
            return &(meshes[i]);
        }

        void Update ()
        {
            generateDistMeshes();
            generateDOFsHandler();
            updated = true;
        }
        Teuchos::RCP<hiperlife::HiPerProblem> generateHiPerProblem(Teuchos::RCP<hiperlife::ParamStructure> paramStr, int gPts=1);

        std::vector<Teuchos::RCP<hiperlife::DistributedMesh>> disMeshes;
        std::vector<Teuchos::RCP<hiperlife::DOFsHandler>> fields;

        void saveSolution(string filename, Teuchos::RCP<hiperlife::ParamStructure> paramStr, double simtime);

    private:

        void generateDistMeshes();
        void generateDOFsHandler();

        bool updated = false;

};

void LS_K (hiperlife::FillStructure& fillStr);
double Getloadstep(double t, double inflation_begin, double inflation_end, double deflation_begin, double deflation_end, double push_mag, double pull_mag);
inline std::string first_numberstring(std::string const & str)
{

    std::size_t const n = str.find_first_of("0123456789");

    if (n != std::string::npos)

    {

        std::size_t const m = str.find_first_not_of("0123456789", n);

        return str.substr(n, m != std::string::npos ? m-n : m);

    }

    return std::string();

}


inline void CalcdGIGRef( double* dGIGRef, double* iTI, double &p_kI, double* Voigt )
{
    for ( int n = 0; n < 3; n++ )
    {
        for (int a = 0; a < 2; a++)
        {
            for (int b = 0; b < 2; b++)
            {
                dGIGRef[4 * n + 2 * a + b] = 0.0;
                for (int c = 0; c < 2; c++)
                {
                    for (int d = 0; d < 2; d++)
                    {
                        dGIGRef[4 * n + 2 * a + b] += iTI[2 * a + c] * iTI[2 * b + d] * p_kI * Voigt[4 * n + 2 * c + d];
                    }
                }
            }
        }
    }
}

inline void JacobianGradient (double *dJ, double *dp_k, double *phi_1, double *phi_2, double *normal)
{
    double cross1[3], cross2[3];

    //cross1 = x_{,1} x n  (3D vector)
    hiperlife::Math::Cross(cross1, phi_1, normal);
    //cross2 = x_{,2} x n (3D vector)
    hiperlife::Math::Cross(cross2, phi_2, normal);


    // dJ = p_{,1}(x_{,2} x n) - p_{,2}(x_{,1} x n)
    hiperlife::Array::Fill3D(dJ, 0.0);
    hiperlife::Math::AXPY(dJ, 3, dp_k[0], cross2);
    hiperlife::Math::AXPY(dJ, 3,-dp_k[1], cross1);

    return;
}

inline void NormalGradient (double *dnormal, double *dp_k, double *phi_1, double *phi_2, double *normal, double J, double *dJ)
{
    double Iphi_1[9], Iphi_2[9];

    hiperlife::Array::Fill(dnormal, 9, 0.0);

    hiperlife::Math::AXPY(&dnormal[0], 3, -normal[0]/J, dJ);
    hiperlife::Math::AXPY(&dnormal[3], 3, -normal[1]/J ,dJ);
    hiperlife::Math::AXPY(&dnormal[6], 3, -normal[2]/J ,dJ);

    Iphi_1[0] =  0.0;
    Iphi_1[1] = -phi_1[2];
    Iphi_1[2] =  phi_1[1];

    Iphi_1[3] =  phi_1[2];
    Iphi_1[4] =  0.0;
    Iphi_1[5] = -phi_1[0];

    Iphi_1[6] = -phi_1[1];
    Iphi_1[7] =  phi_1[0];
    Iphi_1[8] =  0.0;

    Iphi_2[0] =  0.0;
    Iphi_2[1] = -phi_2[2];
    Iphi_2[2] =  phi_2[1];

    Iphi_2[3] =  phi_2[2];
    Iphi_2[4] =  0.0;
    Iphi_2[5] = -phi_2[0];

    Iphi_2[6] = -phi_2[1];
    Iphi_2[7] =  phi_2[0];
    Iphi_2[8] =  0.0;

    hiperlife::Math::AXPY(dnormal, 9, -dp_k[0]/J, Iphi_2);
    hiperlife::Math::AXPY(dnormal, 9,  dp_k[1]/J, Iphi_1);
}

inline void NormalHessian (double *ddnormal, double *dp_kI, double *dp_kJ, double J, double *dJI, double *dJJ, double *ddJ, double *normal, double *dnormalI, double *dnormalJ )
{
    double aux;

    hiperlife::Array::Fill(ddnormal, 27, 0.0);

    aux = -1.0/J * normal[0];
    hiperlife::Math::AXPY(&ddnormal[0], 9, aux, ddJ);
    aux = -1.0/J * normal[1];
    hiperlife::Math::AXPY(&ddnormal[9], 9, aux, ddJ);
    aux = -1.0/J * normal[2];
    hiperlife::Math::AXPY(&ddnormal[18], 9,aux,ddJ);

    aux = -1.0/J * dnormalI[0];
    hiperlife::Math::AXPY(&ddnormal[0], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[1];
    hiperlife::Math::AXPY(&ddnormal[3], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[2];
    hiperlife::Math::AXPY(&ddnormal[6], 3, aux, dJJ);

    aux = -1.0/J * dnormalI[3];
    hiperlife::Math::AXPY(&ddnormal[9], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[4];
    hiperlife::Math::AXPY(&ddnormal[12], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[5];
    hiperlife::Math::AXPY(&ddnormal[15], 3, aux, dJJ);

    aux = -1.0/J * dnormalI[6];
    hiperlife::Math::AXPY(&ddnormal[18], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[7];
    hiperlife::Math::AXPY(&ddnormal[21], 3, aux, dJJ);
    aux = -1.0/J * dnormalI[8];
    hiperlife::Math::AXPY(&ddnormal[24], 3, aux, dJJ);

    aux = -1.0/J * dJI[0];
    hiperlife::Math::AXPY(&ddnormal[0], 3, aux, &dnormalJ[0]);
    aux = -1.0/J * dJI[1];
    hiperlife::Math::AXPY(&ddnormal[3], 3, aux, &dnormalJ[0]);
    aux = -1.0/J * dJI[2];
    hiperlife::Math::AXPY(&ddnormal[6], 3, aux, &dnormalJ[0]);

    aux = -1.0/J * dJI[0];
    hiperlife::Math::AXPY(&ddnormal[9], 3, aux, &dnormalJ[3]);
    aux = -1.0/J * dJI[1];
    hiperlife::Math::AXPY(&ddnormal[12], 3, aux, &dnormalJ[3]);
    aux = -1.0/J * dJI[2];
    hiperlife::Math::AXPY(&ddnormal[15], 3, aux, &dnormalJ[3]);

    aux = -1.0/J * dJI[0];
    hiperlife::Math::AXPY(&ddnormal[18], 3, aux, &dnormalJ[6]);
    aux = -1.0/J * dJI[1];
    hiperlife::Math::AXPY(&ddnormal[21], 3, aux, &dnormalJ[6]);
    aux = -1.0/J * dJI[2];
    hiperlife::Math::AXPY(&ddnormal[24], 3, aux, &dnormalJ[6]);


    aux = 1.0/J * ( dp_kI[0] * dp_kJ[1] - dp_kI[1] * dp_kJ[0] );
    ddnormal[5]  += aux;
    ddnormal[7]  -= aux;
    ddnormal[11] -= aux;
    ddnormal[15] += aux;
    ddnormal[19] += aux;
    ddnormal[21] -= aux;
}

inline void JacobianHessian (double *ddJ, double *dp_kI, double *dp_kJ, double *phi_1, double *phi_2, double *normal, double *dnormalJ)
{
    double aux[3];
    double Inormal[9];
    double cross1[9], cross2[9];

    hiperlife::Array::Fill(ddJ, 9, 0.0);

    Inormal[0] =  0.0;
    Inormal[1] = -normal[2];
    Inormal[2] =  normal[1];

    Inormal[3] =  normal[2];
    Inormal[4] =  0.0;
    Inormal[5] = -normal[0];

    Inormal[6] = -normal[1];
    Inormal[7] =  normal[0];
    Inormal[8] =  0.0;

    hiperlife::Math::AXPY(ddJ, 9, dp_kI[1]*dp_kJ[0]-dp_kI[0]*dp_kJ[1], Inormal);

    aux[0] = dnormalJ[0]; aux[1] = dnormalJ[3]; aux[2] = dnormalJ[6];
    hiperlife::Math::Cross(&cross1[0], phi_1, aux);
    hiperlife::Math::Cross(&cross2[0], phi_2, aux);
    aux[0] = dnormalJ[1]; aux[1] = dnormalJ[4]; aux[2] = dnormalJ[7];
    hiperlife::Math::Cross(&cross1[3], phi_1, aux);
    hiperlife::Math::Cross(&cross2[3], phi_2, aux);
    aux[0] = dnormalJ[2]; aux[1] = dnormalJ[5]; aux[2] = dnormalJ[8];
    hiperlife::Math::Cross(&cross1[6], phi_1, aux);
    hiperlife::Math::Cross(&cross2[6], phi_2, aux);

    hiperlife::Math::AXPY(ddJ, 9,  dp_kI[0], cross2);
    hiperlife::Math::AXPY(ddJ, 9, -dp_kI[1], cross1);

}

// Calculate the product of dN with a vector u: dN.u
inline void NormalGradientU (double *udn, double *dp_k, double *phi_1, double *phi_2, double *normal, double iJac, double *u)
{
    double dJ[3];
    double cross1[3], cross2[3];

    hiperlife::Array::Fill3D(udn, 0.0);

    JacobianGradient(dJ, dp_k, phi_1, phi_2, normal);

    // cross1 = x_{,1} x u
    hiperlife::Math::Cross(cross1, phi_1, u);
    // cross2 = x_{,2} x u
    hiperlife::Math::Cross(cross2, phi_2, u);

    // udn = -iJac * ( u·n·dJ - ( p_{,1} x_{,2} x u - p_{,2} x_{,1} x u) )
    hiperlife::Math::AXPY(udn, 3, -iJac*hiperlife::Math::Dot3D(u, normal), dJ);
    hiperlife::Math::AXPY(udn, 3,  iJac*dp_k[0],                cross2);
    hiperlife::Math::AXPY(udn, 3, -iJac*dp_k[1],                cross1);

}

// Calculate the derivative of the metric tensor with respect to a node coordinates
inline void MetricGradient  (double *dg, double *dp_k, double *phi_1, double *phi_2)
{

    // dg_{ab} = p_{,a} x_{,b} + p_{,b} x_{,a} (12 components, 2x2 tensor and 3 spatial components)
    for (int j = 0; j < 3; j++)
    {
        dg[4*j+0]   = 2.0 * dp_k[0] * phi_1[j];
        double aux  = dp_k[0] * phi_2[j] + dp_k[1] * phi_1[j];
        dg[4*j+1]   = aux;
        dg[4*j+2]   = aux;
        dg[4*j+3]   = 2.0 * dp_k[1] * phi_2[j];
    }

}

// Calculate the derivative of the metric tensor with respect to a node coordinates
inline void MetricHessian(double *ddg, double *dp_kI, double *dp_kJ)
{
    // dg_{ab} = p_{,a} p_{,b} + p_{,b} p_{,a} (should have 3*3*4 components but they are 3 diagonal blocks that are equal + 0's)
    ddg[0]   = 2.0 * dp_kI[0] * dp_kJ[0];
    double aux = dp_kI[0] * dp_kJ[1] + dp_kI[1] * dp_kJ[0];
    ddg[1]   = aux;
    ddg[2]   = aux;
    ddg[3]   = 2.0 * dp_kI[1] * dp_kJ[1];
}



#endif //AUXCORTEXMECHANICS_H
