/*
    *******************************************************************************
    Copyright (c) 2017-2025
    Authors/Contributors: Adam Ouzeri, Sohan Kale, Alejandro Torres-Sánchez, Daniel Santos-Oliván
    *******************************************************************************
    This file is part of agvm_domes
    Project homepage: https://github.com/aouzeri/agvm_domes
    Distributed under the GNU General Public License, see the accompanying
    file LICENSE or https://opensource.org/license/gpl-3-0.
    *******************************************************************************
*/

#include "AuxVXDomes.h"
#include "hl_LocMongeParam.h"
#include "hl_StructMeshGenerator.h"
#include "hl_BasicMeshGenerator.h"
#include "hl_UnstructVtkMeshGenerator.h"
#include "hl_MeshLoader.h"

#include <iostream>
#include <vector>

#include <fstream>
#include <string>
#include <iterator>

using namespace std;
using namespace hiperlife;


void faceMesh::loadMesh(string fileName)
{
    ifstream infile(fileName);
    if(!infile.is_open())
        throw runtime_error("loadFaceMesh: The input file (" + fileName + ") does not exist.");

    string line;

	// Getting the facetype
    getline(infile, line);
    faceType  = stoi(line);

    if ( faceType > 6 or faceType < 0 )
        throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has a incorrent face ID (" + to_string(faceType) + ".");

	// Getting the neighbor faceIDs
    getline(infile, line);

    std::istringstream iss(line);
    vector<string> tokens{istream_iterator<string>{iss},
                          istream_iterator<string>{}};
    for ( auto t: tokens)
        nBorFaces.push_back(stoi(t));
        
        
    // Getting the number of points in the face     
    getline(infile, line);
    nPts  = stoi(line);
    getline(infile, line);
    nElem = stoi(line);

	// Getting the nodal coordinates
    x_nodes.resize(nPts*nDim);
    cn_nodes.resize(nPts*nDim);

    for ( int i = 0; i < nPts; i++ )
    {
        getline(infile, line);

        if ( infile )
        {
            if (line.empty())
                throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has a blank line.");

            std::istringstream iss(line);
            vector<string> tokens{istream_iterator<string>{iss},
                                  istream_iterator<string>{}};

            if (tokens.size() != nDim)
                throw runtime_error("loadFaceMesh: The input file (" + fileName +
                                    ") has a line with an inconsistent number of nodal coordinates.");

            for ( int n = 0; n < nDim; n++)
                x_nodes[nPts*n+i] = stod(tokens[n]);

            // Calculating the centroid position at that mesh
            centroidPositionInFaceMesh[0] += stod(tokens[0])/nPts;
            centroidPositionInFaceMesh[1] += stod(tokens[1])/nPts;
            centroidPositionInFaceMesh[2] += stod(tokens[2])/nPts;

        }
        else
            throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has finished before the end of reading. Inconsistent data has been provided.");
    }

	// Getting the connectivity matrix
    for ( int e = 0; e < nElem; e++ )
    {
        getline(infile, line);

        if ( infile )
        {
            if (line.empty())
                throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has a blank line.");

            std::istringstream iss(line);
            vector<string> tokens{istream_iterator<string>{iss},
                                  istream_iterator<string>{}};

            if (tokens.size() != nVert)
                throw runtime_error("loadFaceMesh: The input file (" + fileName +
                                    ") has a line with an inconsistent number of nodes per element.");

            for ( int n = 0; n < nVert; n++)
                connec.push_back(stoi(tokens[n]));
        }
        else
            throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has finished before the end of reading. Inconsistent data has been provided.");
    }

	// Getting the master-slave relationship: faceID of master, (local) nodeID of master, additive component to relationship (0 if they should be equal)
    masterNodes.resize(nPts*nDim); // Storing separately for each DOF
    masterNodesPerFlag.resize(nPts*nDim); // For storing the periodic BC flags
    for (int i = 0; i < nPts; i++)
    {
        for (int dof = 0; dof < nDim; dof++)
        {
            getline(infile, line);

            if (infile)
            {
                if (line.empty())
                    continue;

                std::istringstream iss(line);
                vector<string> tokens{istream_iterator<string>{iss},
                                      istream_iterator<string>{}};

                if (tokens.size() % 3 != 0)
                {
                    cout << fileName << ":"<< tokens.size() << "," << line.empty() << ", [" << line << "]" << endl;
		            throw runtime_error("loadFaceMesh: The input file (" + fileName +
                                      ") has a line with an inconsistent number of inputs for LC.");
                }

                for (int n = 0; n < tokens.size() / 3; n++)
                {
                    masterNodes[nPts*dof + i].push_back({stoi(tokens[3 * n]), stoi(tokens[3 * n + 1])});
                    masterNodesPerFlag[nPts*dof + i].push_back(stod(tokens[3 * n + 2]));
                }
            }
            else
            {
                throw runtime_error("loadFaceMesh: The input file (" + fileName +
                                    ") has finished before the end of reading. Inconsistent data has been provided.");
            }
        }
    }

	// Getting constraints
    for ( int i = 0; i < nPts; i++ )
    {
        getline(infile, line);

        if ( infile )
        {
            if (line.empty())
                throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has a blank line.");

            std::istringstream iss(line);
            vector<string> tokens{istream_iterator<string>{iss},
                                  istream_iterator<string>{}};

            if (tokens.size() != nDim)
                throw runtime_error("loadFaceMesh: The input file (" + fileName +
                                    ") has a line with an inconsistent number of constraints.");

            for ( int n = 0; n < nDim; n++)
                cn_nodes[nPts*n+i] = bool(stoi(tokens[n]));
        }
        else
            throw runtime_error("loadFaceMesh: The input file (" + fileName + ") has finished before the end of reading. Inconsistent data has been provided.");
    }


    idN_creases.resize(nPts);
    // tag on the node for boundary (mine is always)
    fill(idN_creases.begin(),idN_creases.end(),0);
    

}

void tissueMesh::loadMesh( string fileName)
{
    string infoFileName = fileName + "_info.txt";

    ifstream inInfoFile(infoFileName);

    string line;
    getline(inInfoFile, line);
    nFaces = stoi(line);
    nPointsinFace.resize(nFaces,0);

    setNItem(nFaces);
    DistributedData::Update();

    vector<int > loc_nPointsinFace(nFaces,0);


    for ( int n = 0; n < nFaces; n++ )
    {
        if ( isItemInPart(n) )
        {
            faceMesh fmesh;
            fmesh.loadMesh(fileName + "_" + to_string(n)+".txt");
            loc_nPointsinFace[n] = fmesh.nPts;
            meshes.push_back(fmesh);
        }
    }
    
    MPI_Allreduce(loc_nPointsinFace.data(),nPointsinFace.data(),nFaces,MPI_INT,MPI_MAX,_comm);

}

void tissueMesh::generateDistMeshes() {
    if(myRank() == 0)
        cout << "Generating distributed mesh..." << endl;
    
    for (int i = 0; i < nItem(); i++) {


        SmartPtr<BasicMeshGenerator> _mesh = Create<BasicMeshGenerator>();
        _mesh->setBasisFuncType(BasisFuncType::Linear);
        _mesh->setBasisFuncOrder(1);
        _mesh->setMeshType(MeshType::Parallel);
        _mesh->setElemType(ElemType::Triang);

        if (isItemInPart(i)) {
            faceMesh *fmesh = &meshes[locIdx(i)];
            _mesh->setLocNPts(fmesh->nPts);
            _mesh->setLocNElem(fmesh->nElem);
            _mesh->setDistCoords(fmesh->x_nodes);
            _mesh->setDistConnec(fmesh->connec);
            _mesh->setDistCreases(fmesh->idN_creases);
            _mesh->_eflags.resize(fmesh->nElem);
        }

        //Distribute the mesh
        SmartPtr<DistributedMesh> _disMesh = Create<DistributedMesh>();
        _disMesh->setMesh(_mesh);
        _disMesh->setBalanceMesh(false);
        _disMesh->setElementLocatorEngine(ElementLocatorEngine::None);
        _disMesh->Update();

        //Save it in vector
        disMeshes.push_back(_disMesh);
    }

    if(myRank() == 0)
        cout << "Finished generating distributed mesh..." << endl;
}

void tissueMesh::generateDOFsHandler() // Now called DOFsHandler
{
        if(myRank() == 0)
        cout << "Generating fields..." << endl;

    int  DOF = meshes[0].DOF;
    int gDOF = meshes[0].gDOF;
    int gAux = meshes[0].gAux;

    for ( int i = 0; i < nItem(); i++)
    {
		// Node DOFs will depend on the structured mesh
        SmartPtr<DOFsHandler> _fields = Create<DOFsHandler>(disMeshes[i]);

        if ( isItemInPart(i))
        {
            int locI = locIdx(i);

            switch ( meshes[locI].faceType )
            {
                case 0:
                    _fields->setNameTag("apical"+to_string(i));
                    break;
                case 1:
                    _fields->setNameTag("basal"+to_string(i));
                    break;
                case 2:
                    _fields->setNameTag("lateral"+to_string(i));
                    break;
                case 3:
                    _fields->setNameTag("lL"+to_string(i));
                    break;
                case 4:
                    _fields->setNameTag("lR"+to_string(i));
                    break;
                case 5:
                    _fields->setNameTag("lB"+to_string(i));
                    break;
                case 6:
                    _fields->setNameTag("lT"+to_string(i));
                    break;
            }
        }
        _fields->setNumDOFs(DOF);
        _fields->setNumNodeAuxF(3);   // to store the nodal forces
        _fields->Update();

		// gDOFs will depend on the global constraint
        SmartPtr<DistributedMesh> _gloMesh  = Create<DistributedMesh>();
        _gloMesh->setMeshRelation(MeshRelation::GlobConstr, disMeshes[i]);
        _gloMesh->setMasterProcessor(getItemPartition(i));
        _gloMesh->setBalanceMesh(false);
        _gloMesh->setElementLocatorEngine(ElementLocatorEngine::None);
        _gloMesh->Update();

        SmartPtr<DOFsHandler> _gfields = Create<DOFsHandler>(_gloMesh);
        _gfields->setNameTag("gcons"+to_string(i));
        _gfields->setNumDOFs(gDOF);
        _gfields->setNumNodeAuxF(gAux); // to store the local-monge basis (2 basis vector 3D components)
        _gfields->Update();

        if(isItemInPart(i))
        {
            int locI      = locIdx(i);

            meshes[locI].nodeDOFs =  _fields->nodeDOFs;
            meshes[locI].auxF     =  _fields->nodeAuxF;
            meshes[locI].gDOFS    = _gfields->nodeDOFs;
            meshes[locI].auxV     = _gfields->nodeAuxF;

            // Read from input files to do a restart
            if ( meshes[locI].u_nodes.size() > 0 )
            {
                for ( int j = 0; j < meshes[locI].nPts; j++ )
                {
                    for ( int n = 0; n < meshes[locI].DOF; n++ )
                    {
                        meshes[locI].nodeDOFs->setValue(n,j,IndexType::Local,meshes[locI].u_nodes[j*DOF+n]);
                    }
                }
            }
            if ( meshes[locI].globalU.size() > 0 )
            {
                for ( int n = 0; n < meshes[locI].gDOF; n++ )
                {
                    meshes[locI].gDOFS->setValue(n,0,IndexType::Local,meshes[locI].globalU[n]);
                }
            }
            if ( meshes[locI].globalV.size() > 0 )
            {
                for ( int n = 0; n < meshes[locI].gAux; n++ )
                {
                    meshes[locI].auxV->setValue(n,0,IndexType::Local,meshes[locI].globalV[n]);
                }
            }

            int loc_nPts = meshes[locI].nPts;

            // Initialize the DOFs and local-Monge basis to the default values //FIXME: MOVE TO FUNCTION ?
                //Calculate basis
                int n0e = disMeshes[i]->_connec->getNbor(0, 0, IndexType::Local);
                int n1e = disMeshes[i]->_connec->getNbor(0, 1, IndexType::Local);
                int n2e = disMeshes[i]->_connec->getNbor(0, 2, IndexType::Local);


                vector<double> x0e = disMeshes[i]->nodeCoords(n0e);
                vector<double> x1e = disMeshes[i]->nodeCoords(n1e);
                vector<double> x2e = disMeshes[i]->nodeCoords(n2e);


                double g1e[3] = {};
                double g2e[3] = {};

                Array::Copy(g1e, 3, x1e.data());
                Math::AXPY(g1e, 3, -1.0, x0e.data());
                Array::Copy(g2e, 3, x2e.data());
                Math::AXPY(g2e, 3, -1.0, x0e.data());

                Math::AX(g1e, 3, 1.0 / Math::Norm3D(g1e));

                Math::AXPY(g2e, 3, -Math::Dot3D(g2e, g1e), g1e);
                Math::AX(g2e, 3, 1.0 / Math::Norm3D(g2e));

                for (int n = 0; n < 3; n++) {
                    _gfields->nodeAuxF->setValue(0 + n, 0, IndexType::Global, g1e[n]);
                    _gfields->nodeAuxF->setValue(3 + n, 0, IndexType::Global, g2e[n]);
                }


                _gfields->nodeDOFs->setValue(0, 0, IndexType::Local, 1.0); //G11
                _gfields->nodeDOFs->setValue(1, 0, IndexType::Local, 0.0); //G12
                _gfields->nodeDOFs->setValue(2, 0, IndexType::Local, 1.0); //G22

                _gfields->nodeDOFs->setValue(3, 0, IndexType::Local, 0.0);  // Cell pressure
                _gfields->nodeDOFs->setValue(4, 0, IndexType::Local, 0.03); // Density
                _gfields->nodeDOFs->setValue(5, 0, IndexType::Local, 0.0); // Lumen pressure

                // Initialize nodeDOFs and nodeAUX values
                for (int j = 0; j < loc_nPts; j++) {
                    vector<double> x = _fields->mesh->nodeCoords(j, IndexType::Global);
                    _fields->nodeDOFs->setValue(0, j, IndexType::Global, x[0]);
                    _fields->nodeDOFs->setValue(1, j, IndexType::Global, x[1]);
                    _fields->nodeDOFs->setValue(2, j, IndexType::Global, x[2]);

                    // to store nodal forces
                    _fields->nodeAuxF->setValue(0, j, IndexType::Global, 1.0);
                    _fields->nodeAuxF->setValue(1, j, IndexType::Global, 1.0);
                    _fields->nodeAuxF->setValue(2, j, IndexType::Global, 1.0);
                }

            // Set dirichlet constraints
            for ( int j = 0; j < loc_nPts; j++ )
            {
                if (setConstraintflag == 1)
                {
                    for ( int n = 0; n < meshes[locI].nDim; n++ )
                    {
                        if ( meshes[locI].cn_nodes[meshes[locI].nPts*n+j] > 0)
                        {
                            _fields->setConstraint(n, j, IndexType::Local, 0.0); // n = DOF, j = local node index (of that face/mesh i)
                        }
                    }
                }
            }
        }

        _fields->nodeDOFs0->setValue(_fields->nodeDOFs);
        _gfields->nodeDOFs0->setValue(_gfields->nodeDOFs);

        _fields->UpdateGhosts();
        _gfields->UpdateGhosts();

        fields.push_back(_fields); // for nodes DOFs (nodes position)
        fields.push_back(_gfields); // for global DOFs (density, pressure, etc.)
    }
}

void tissueMesh::saveSolution(string filename, SmartPtr<ParamStructure> paramStr, double simtime)
{
    std::vector<int>  offNodes;
    std::vector<int> goffNodes;

    offNodes.resize(loc_nItem()+1);
    goffNodes.resize(loc_nItem()+1); // global offset

    int loc_nPts{};
    int loc_nElem{};

    int n = 0;
    offNodes[0] = 0;
    for ( int i = 0; i < loc_nItem(); i++ )
        offNodes[i+1] = offNodes[i] + meshes[i].nPts;

    loc_nPts = offNodes[loc_nItem()];

    DistributedData partoffNodes;
    partoffNodes.setLocNItem(loc_nPts);
    partoffNodes.Update();

    for ( int i = 0; i < loc_nItem()+1; i++ )
        goffNodes[i] = offNodes[i] + partoffNodes.offs_nItem();

    vector<int> sconnec={};
    vector<double> sx_nodes={};
    vector<int> screases={};
    vector<double> su_nodes={};


    sx_nodes.resize(loc_nPts*meshes[0].nDim);

    for ( int i = 0; i < loc_nItem(); i++ )
    {
        auto _mesh        = meshes[i];
        auto _connec      = _mesh.connec;
        auto _x_nodes     = _mesh.x_nodes;
        auto _idN_creases = _mesh.idN_creases;
        int  _nVert       = _mesh.nVert;
        int  _nDim        = _mesh.nDim;
        int  _nPts        = _mesh.nPts;
        int  _nElem       = _mesh.nElem;

        loc_nElem += _nElem;
        for ( int e = 0; e < _nElem; e++ )
        {
            for ( int n = 0; n < _nVert; n++)
                sconnec.push_back(_connec[e*_nVert+n]+goffNodes[i]);
        }

        for ( int j = 0; j <_nPts; j++ )
        {
            for ( int n = 0; n < _nDim; n++ )
                sx_nodes[loc_nPts*n+j+offNodes[i]] = _x_nodes[_nPts*n+j];

        }

        screases.insert(screases.end(), _idN_creases.begin(), _idN_creases.end());
    }

    SmartPtr<StructMeshGenerator> mesh = Create<StructMeshGenerator>();
    mesh->setMesh(ElemType::Triang,BasisFuncType::Linear,1);

    mesh->setNDim(meshes[0].nDim);
    mesh->_loc_nPts  = loc_nPts;
    mesh->_loc_nElem = loc_nElem;

    mesh->_connec   = sconnec;
    mesh->_x_nodes  = sx_nodes;
    mesh->_creases  = screases;
    mesh->_eflags.resize(loc_nElem);

    int sum_nElem;
    Barrier();
    MPI_Allreduce(&loc_nElem, &sum_nElem, 1, MPI_INT, MPI_SUM, _comm);

    mesh->_nPts  = 1;
    mesh->_nElem = sum_nElem;
    mesh->_nVert = 1;
    mesh->_pDim  = 2;

    SmartPtr<DistributedMesh> disMesh = Create<DistributedMesh>();

    disMesh->setMesh(mesh);
    disMesh->setBalanceMesh(false);
    disMesh->setElementLocatorEngine(ElementLocatorEngine::None);
    disMesh->Update();

    // Generate dummy DOFsHandler
    SmartPtr<DOFsHandler> dofsHand = Create<DOFsHandler>(disMesh);
    dofsHand->setNumDOFs(meshes[0].DOF);
    dofsHand->setNumNodeAuxF(2); 
    dofsHand->Update();

    for ( int i = 0; i < loc_nItem(); i++ )
    {
        auto _mesh        = meshes[i];
        int  _nPts        = _mesh.nPts;
        int  _DOF         = _mesh.DOF;
	
        for ( int j = 0; j <_nPts; j++ )
        {
			// writing the nodes deformed coordinates in the dummy DOFsHandler
            for ( int n = 0; n < _DOF; n++ )
            {
                double val = meshes[i].nodeDOFs->getValue(n,j,IndexType::Local);
                dofsHand->nodeDOFs->setValue(n,j+goffNodes[i],IndexType::Global,val);
            }
            
            // writing density in the dummy DOFsHandler
            double density = meshes[i].gDOFS->getValue(4, 0, IndexType::Local);
            dofsHand->nodeAuxF->setValue(0, j + goffNodes[i], IndexType::Global,density);

            // saving face type (0 = apical, 1 = basal, > 1 = lateral )
            dofsHand->nodeAuxF->setValue(1, j + goffNodes[i], IndexType::Global, meshes[i].faceType); 
        }
    
    }

    // Printing the .vtu and .vtm files
    dofsHand->printFileVtk(filename,true, simtime);
}

SmartPtr<HiPerProblem> tissueMesh::generateHiPerProblem(SmartPtr<ParamStructure> paramStr, int gPts)
{
    if ( !updated )
        throw runtime_error("tissueMesh has not been previously updated. Can't generate HiPerProblem out of nothing!");

    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>("HiperObject");
    hiperProbl->setParameterStructure(paramStr);
    hiperProbl->setDOFsHandlers(fields);
    hiperProbl->setConsistencyCheckType(ConsistencyCheckType::None);
    
    for ( int i = 0; i < nItem(); i++)
    {
        // Doing an integration over each face (so defining a new integration for each face)
        string integName = "Integ"+to_string(i);

        hiperProbl->setIntegration(integName, {fields[2*i+0]->nameTag(),fields[2*i+1]->nameTag()});
        hiperProbl->setCubatureGauss(integName, gPts);
        hiperProbl->setElementFillings(integName, LS_K);
    }

    int loc_nFaces = loc_nItem();
    
    // Setting linear constraints (master-slave)
    if (myRank() == 0)
	    cout << "Setting the DOF linear constraints" << endl;

    for (int n = 0; n < loc_nFaces; n++)
    {
        int loc_nPts = meshes[n].nPts;

        for (int i = 0; i < loc_nPts; i++)
        {
            int nMasters = 0;
            int nDim     = meshes[n].nDim;

            for (int dof = 0; dof < nDim; dof++)
            {
                nMasters = (meshes[n].masterNodes[loc_nPts*dof + i].size());

                if (nMasters > 0)
                {
                    double weight = 1.0 / nMasters;
                    for (int masterID = 0; masterID < nMasters; masterID++)
                    {
                        auto master = meshes[n].masterNodes[loc_nPts*dof + i][masterID];
                        int face = master.first;
                        int node = master.second;
                        // Setting linear constraints between different faces
                        hiperProbl->setLinearConstraint({2 * (n + offs_nItem()), dof, i, IndexType::Global}, 
                                                    {2 * face, dof, node, IndexType::Global}, weight, 0.0);
                    }
                }
            }
        }

        // Set cell pressures to be the same for all faces belonging to the same cell
        int maxNBorface = *max_element(meshes[n].nBorFaces.begin(), meshes[n].nBorFaces.end());
        int face = n + offs_nItem();
        if (face < maxNBorface)
        {
              hiperProbl->setLinearConstraint({2 * (n + offs_nItem()) + 1, 3, 0, IndexType::Global},
                                           {2 * maxNBorface + 1, 3, 0, IndexType::Global}, 1.0, 0.0);
        }

        // Set lumen pressure to be the same as defined in the gDOFs of the 0th face
        if (face > 0)
         {
            hiperProbl->setLinearConstraint({2 * (n + offs_nItem()) + 1, 5, 0, IndexType::Global},
                                           {1, 5, 0, IndexType::Global}, 1.0, 0.0);
         }

    }


    if (myRank() == 0)
	    cout << "Updating hiperproblem..." << endl;

    hiperProbl->Update();

    if (myRank() == 0)
	    cout << "Updating ghosts..." << endl;

    hiperProbl->UpdateGhosts();

    return hiperProbl;
}


void tissueMesh::setNumberCells(const int numberCells)
{
    if(numberCells < 1)
        throw runtime_error("There are only " + to_string(numberCells) + " cells in the tissue.");

    nCells = numberCells;
    NeighboursCells.resize(nCells);
    facesInCell.resize(nCells);
}

//FIXME: parallelise ?
void tissueMesh::readNeighbouringCells(const string fileName)
{
    ifstream infile(fileName);
    if(!infile.is_open())
        throw runtime_error("readNeighbouringCellsfromFile: The input file (" + fileName + ") does not exist.");

    string line;

    // Getting the neighbor cellIDs
    for(int i = 0; i < nCells; i++)
    {
        getline(infile, line);

        std::istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};
        for ( auto t: tokens)
            NeighboursCells[i].push_back(stoi(t));

    }

}

void tissueMesh::readfaceIDsinCells(const string fileName)
{
    ifstream infile(fileName);
    if(!infile.is_open())
        throw runtime_error("readfaceIDsingCells: The input file (" + fileName + ") does not exist.");

    string line;

    // Getting the faceIDs in cellIDs
    for(int c = 0; c < nCells; c++)
    {
        getline(infile, line);

        std::istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};
        for ( auto t: tokens)
            facesInCell[c].push_back(stoi(t));

    }
}

void tissueMesh::setNPointsInCell(SmartPtr<ParamStructure> paramStr)
{
 // Assumes all lateralfaces have same number of points

	vector <int> nPointsinCell(nCells,0);
	vector <int> nPointsinCellAtFace(nFaces,0);
	for(int c = 0; c < nCells; c++) {
        	int nFacesinCell = facesInCell[c].size();
		bool countedLateralFace = false;
		int nVertEdgepoints = 0;
		int nHorEdgepoints = 0;
        	for (int m = 0; m < nFacesinCell; m++) {
            		int globfaceID = facesInCell[c][m];
            		
			if (isItemInPart(globfaceID)) {
	
                		int locI = locIdx(globfaceID );
				int loc_nPts = meshes[locI].nPts;
				int faceType = meshes[locI].faceType;
				if(faceType < 2)
					{nPointsinCell[c] += loc_nPts;}
				else if (!countedLateralFace)
				{
					countedLateralFace = true;
					double xlocation0 = meshes[locI].x_nodes[loc_nPts*0 + 0];
					double ylocation0 = meshes[locI].x_nodes[loc_nPts*1 + 0];
					double zlocation0 = meshes[locI].x_nodes[loc_nPts*2 + 0];

					for(int i = 0; i < loc_nPts; i++)
					{
						if(abs(meshes[locI].x_nodes[loc_nPts*0 + i] - xlocation0) < 1e-10 and abs(meshes[locI].x_nodes[loc_nPts*1 + i] - ylocation0) < 1e-10)
							{nVertEdgepoints += 1;}
						
						if(abs(meshes[locI].x_nodes[loc_nPts*2 + i] - zlocation0) < 1e-10)
							{nHorEdgepoints += 1;}
	
					}
				}
            						}
        						}
		nPointsinCell[c] += (nFacesinCell - 2)*(nVertEdgepoints - 2)*(nHorEdgepoints - 1);
                                      
    					}

	for(int c = 0; c < nCells; c++) 
		for(int f = 0; f < facesInCell[c].size(); f++)
			nPointsinCellAtFace[facesInCell[c][f]] = nPointsinCell[c];

// Sharing the information with all the partition
    MPI_Allreduce(nPointsinCellAtFace.data(), paramStr->i_aux.data(), nItem(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);

}

void tissueMesh::updateFacesCentroids(SmartPtr<ParamStructure> paramStr)
{

    int loc_nFaces = loc_nItem();

    for (int locI = 0; locI < loc_nFaces; locI++) {
        //            cout << "Face: " << face << " is in partition: " << myRank() << endl;

        // Resetting to 0 for new calculations
            meshes[locI].centroidPositionInFaceMesh[0] = 0;
            meshes[locI].centroidPositionInFaceMesh[1] = 0;
            meshes[locI].centroidPositionInFaceMesh[2] = 0;
            double nPoints =  meshes[locI].nPts;

            for (int p = 0; p < nPoints; p++) {
            meshes[locI].centroidPositionInFaceMesh[0] += meshes[locI].nodeDOFs->getValue(0, p, IndexType::Local)/nPoints;
            meshes[locI].centroidPositionInFaceMesh[1] += meshes[locI].nodeDOFs->getValue(1, p, IndexType::Local)/nPoints;
            meshes[locI].centroidPositionInFaceMesh[2] += meshes[locI].nodeDOFs->getValue(2, p, IndexType::Local)/nPoints;

        }
    }
}


void tissueMesh::calculateDistanceBetweenCellCentroids(SmartPtr<ParamStructure> paramStr, SmartPtr<HiPerProblem> hiperProbl)
{

    vector<double > centroidPositionInCell(3*nCells, 0.0);
    vector<double>  CentroidOfCells(nCells*3);

    //Calculating centroid position of each cell
    for(int c = 0; c < nCells; c++) {
        int nFacesinCell = facesInCell[c].size();

        for (int m = 0; m < nFacesinCell; m++) {
            int face = facesInCell[c][m];
            if (isItemInPart(face)) {

                int locI = locIdx(face);
		double nPoints =  meshes[locI].nPts;

                centroidPositionInCell[3*c + 0] += meshes[locI].centroidPositionInFaceMesh[0]/ nFacesinCell;
                centroidPositionInCell[3*c + 1] += meshes[locI].centroidPositionInFaceMesh[1]/ nFacesinCell;
                centroidPositionInCell[3*c + 2] += meshes[locI].centroidPositionInFaceMesh[2]/ nFacesinCell;
            }
        }
    }

    // Sharing the centroidPositioninCell with all the partitions
    MPI_Allreduce(centroidPositionInCell.data(),CentroidOfCells.data(), nCells*3, MPI_DOUBLE, MPI_SUM, _comm);


    // Calculating distances between cell centroids
    for(int c1 = 0; c1 < nCells; c1++) {
        for (int c2 = c1 + 1; c2 < nCells; c2++) {
            double x2 = (CentroidOfCells[3 * c1 + 0] - CentroidOfCells[3 * c2 + 0]) *
                        (CentroidOfCells[3 * c1 + 0] - CentroidOfCells[3 * c2 + 0]);
            double y2 = (CentroidOfCells[3 * c1 + 1] - CentroidOfCells[3 * c2 + 1]) *
                        (CentroidOfCells[3 * c1 + 1] - CentroidOfCells[3 * c2 + 1]);
            double z2 = (CentroidOfCells[3 * c1 + 2] - CentroidOfCells[3 * c2 + 2]) *
                        (CentroidOfCells[3 * c1 + 2] - CentroidOfCells[3 * c2 + 2]);
            double distCells = sqrt(x2 + y2 + z2);


            // Once two cells have been in contact, we check which faces have touched (the ones where the distance between their centroids is the smallest).
            // Then we apply linear constraints on theses two faces to simulate sticking between cells.
            double distCellsThreshold = 1.8;
            if (find(NeighboursCells[c1].begin(), NeighboursCells[c1].end(), c2) == NeighboursCells[c1].end())  // if not a neighbour
                if (distCells >= distCellsThreshold) {
                        // DO NOTHING. SETTING TO FALSE IS DONE IN THE MAIN LOOP
                    } else {


                        // retrieving the centroids of each face in cell1 and cell2
                        vector<double> posList1(facesInCell[c1].size() * 3);
                        vector<double> posList2(facesInCell[c2].size() * 3);

                        for (int f = 0; f < facesInCell[c1].size(); f++) {
                            int face = facesInCell[c1][f];
                            if (isItemInPart(face)) {
                                int locI = locIdx(face);
                                posList1[3 * f + 0] = meshes[locI].centroidPositionInFaceMesh[0];
                                posList1[3 * f + 1] = meshes[locI].centroidPositionInFaceMesh[1];
                                posList1[3 * f + 2] = meshes[locI].centroidPositionInFaceMesh[2];
                            }


                        }

                        for (int f = 0; f < facesInCell[c2].size(); f++) {
                            int face = facesInCell[c2][f];
                            if (isItemInPart(face)) {
                                int locI = locIdx(face);
                                posList2[3 * f + 0] = meshes[locI].centroidPositionInFaceMesh[0];
                                posList2[3 * f + 1] = meshes[locI].centroidPositionInFaceMesh[1];
                                posList2[3 * f + 2] = meshes[locI].centroidPositionInFaceMesh[2];
                            }
                        }

                        vector<double> CentroidsOfFacesinC1(facesInCell[c1].size() * 3);
                        vector<double> CentroidsOfFacesinC2(facesInCell[c2].size() * 3);

                        // Sharing with all the partitions
                        MPI_Allreduce(posList1.data(), CentroidsOfFacesinC1.data(), facesInCell[c1].size() * 3,
                                      MPI_DOUBLE, MPI_SUM, _comm);
                        MPI_Allreduce(posList2.data(), CentroidsOfFacesinC2.data(), facesInCell[c2].size() * 3,
                                      MPI_DOUBLE, MPI_SUM, _comm);

                        // Calculating all the distances between each face of a cell with the barycenter of the other cell and taking the closest.
			            int face1ID = -1;
                        int face2ID = -1;
                        int locface1ID = -1;
                        int locface2ID = -1;
			            int faceNumberinC1 = -1;
			            int faceNumberinC2 = -1;
                        double mindistbetweenFaceIandFaceII = 1000.0; // high value for comparison


                    double nPointsinCell1 = 0;
                    double nPointsinCell2 = 0;

                    // comparing distance of faces in cell 1 and faces in cell 2
                    for (int faceI = 0; faceI < facesInCell[c1].size(); faceI++) {
                        nPointsinCell1 = nPointsinCell1 + nPointsinFace[facesInCell[c1][faceI]]; //!Adam : counts duplicates

                        for (int faceII = 0; faceII < facesInCell[c2].size(); faceII++) {
                            if(faceI == 0)
                                nPointsinCell2 = nPointsinCell2 + nPointsinFace[facesInCell[c2][faceII]];

                            double x2 = (CentroidsOfFacesinC1[3 * faceI + 0] - CentroidsOfFacesinC2[3 * faceII + 0]) *
                                        (CentroidsOfFacesinC1[3 * faceI + 0] - CentroidsOfFacesinC2[3 * faceII + 0]);
                            double y2 = (CentroidsOfFacesinC1[3 * faceI + 1] - CentroidsOfFacesinC2[3 * faceII + 1]) *
                                        (CentroidsOfFacesinC1[3 * faceI + 1] - CentroidsOfFacesinC2[3 * faceII + 1]);
                            double z2 = (CentroidsOfFacesinC1[3 * faceI + 2] - CentroidsOfFacesinC2[3 * faceII + 2]) *
                                        (CentroidsOfFacesinC1[3 * faceI + 2] - CentroidsOfFacesinC2[3 * faceII + 2]);
                            double distbetweenFaceIandFaceII = sqrt(x2 + y2 + z2);

                            if (distbetweenFaceIandFaceII < mindistbetweenFaceIandFaceII - 1e-6) {
                                face1ID = facesInCell[c1][faceI];
                                face2ID = facesInCell[c2][faceII];
                                locface1ID = faceI;
                                locface2ID = faceII;
                                mindistbetweenFaceIandFaceII = distbetweenFaceIandFaceII;
                            }
                        }
                    }

                    if (abs(mindistbetweenFaceIandFaceII) < 1.0){

                        for (int faceI = 0; faceI < facesInCell[c1].size(); faceI++) {
                            int globFace1ID = facesInCell[c1][faceI];
                            paramStr->bparam[globFace1ID] = true;
                            paramStr->x_aux[globFace1ID] = (CentroidOfCells[3 * c2 + 0] - CentroidOfCells[3 * c1 + 0]); // / distCells / distCells / distCells / nPointsinCell1;
                            paramStr->y_aux[globFace1ID] = (CentroidOfCells[3 * c2 + 1] - CentroidOfCells[3 * c1 + 1]); // / distCells / distCells / distCells / nPointsinCell1;
                            paramStr->z_aux[globFace1ID] = (CentroidOfCells[3 * c2 + 2] - CentroidOfCells[3 * c1 + 2]); // / distCells / distCells / distCells / nPointsinCell1;
                        }

                        for (int faceII = 0; faceII < facesInCell[c2].size(); faceII++) {
                            int globFace2ID = facesInCell[c2][faceII];
                            paramStr->bparam[globFace2ID] = true;
                            paramStr->x_aux[globFace2ID] = (CentroidOfCells[3 * c1 + 0] - CentroidOfCells[3 * c2 + 0]);// / distCells / distCells / distCells / nPointsinCell2;
                            paramStr->y_aux[globFace2ID] = (CentroidOfCells[3 * c1 + 1] - CentroidOfCells[3 * c2 + 1]);// / distCells / distCells / distCells / nPointsinCell2;
                            paramStr->z_aux[globFace2ID] = (CentroidOfCells[3 * c1 + 2] - CentroidOfCells[3 * c2 + 2]);// / distCells / distCells / distCells / nPointsinCell2;
                        }

                   }

                }
        }

    }

}


void LS_K (hiperlife::FillStructure& fillStr)
{
    //------------------------------------------------------------------
    // [1] INPUT definition

    //[1.1] Load main subfill structure:
    SubFillStructure& subFill = fillStr[0];

    int nDim  = subFill.nDim; // nDim = 3
    int DOF   = subFill.numDOFs; // DOF = 3                   //Number of degrees of freedom (3)
    int eNN   = subFill.eNN; // eNN = 3                       //Number of neighbors to the element (3 Linear Elements)
    int pDim  = subFill.pDim; // pDim = 2                      //Parametric dimension (2)
    double *xe_nodes  = subFill.nborCoords.data();  //Reference surface positions
    double *ue_nodes  = subFill.nborDOFs.data();    //Degrees of freedom
    double *ue_nodes0 = subFill.nborDOFs0.data();   //Degrees of freedom of previous time-step

    int N_k = DOF*eNN; // N_k = 9 //Row size
    double *fe_nodes  = subFill.nborAuxF.data();  // Auxillary data defined on the nodes

    //FIXME: reverse map?
    std::map<std::string,int>::const_iterator it;
    string key = {};
    for (it = fillStr._dhandsTagIdx.begin(); it != fillStr._dhandsTagIdx.end(); ++it)
    {
        if (it->second == 0) // 0 - points to the master mechanics, 1 - points to the global mechanics
        {
            key = it->first;
            break;
        }
    }

    // Extract the faceID from the key string
    int faceID = std::atoi(first_numberstring(key).c_str());

    //[1.2] Load global DOF subfill structure
    SubFillStructure& g_subFill = fillStr[1];
    int gDOF   = g_subFill.numDOFs; 				//Number of global degrees of freedom 
    int geNN   = g_subFill.eNN;
    int gN_k   = geNN * gDOF;
    double *g_ue_nodes  = g_subFill.nborDOFs.data();
    double *g_ue_nodes0 = g_subFill.nborDOFs0.data();
    double *g_ae_nodes  = g_subFill.nborAuxF.data();

    double pressure = g_ue_nodes[3];
    double dens = g_ue_nodes[4];
    double dens0 = g_ue_nodes0[4];
    double lumpressure = g_ue_nodes[5];

    double conc = fillStr.paramStr->dparam[11];

    //[1.3] Load basis functions and derivatives
    double *p_k  = subFill.nborBasisFunctionsDerivatives(0);
    double *dp_k = subFill.nborBasisFunctionsDerivatives(1);

    //[1.4] Load model parameters
    vector<double>   dparam = fillStr.paramStr->dparam;
    double deltat    = dparam[0];
    double kp        = dparam[1] * deltat;
    double kd        = dparam[2] * deltat;
    double bVisc     = dparam[3];
    double lame1     = dparam[4] * deltat;
    double lame2     = dparam[5] * deltat;
    double eta_f     = dparam[9];
    double decreasecoeff = dparam[12];
    double contact_force_coefficient = dparam[13];

    double vdome         = dparam[10];
    double activ{};
    int    facetag   = 0;

    switch(key[0])
    {
        case 'a':
            activ = dparam[6] * deltat;
            facetag = 0;
            break;
        case 'b':
            activ = dparam[7] * deltat;
            facetag = 1;
            break;
        case 'l':
        // To keep r = gamma_l/gamma_ab = mu_l / mu_ab = eta_l / eta_ab (dparam[8])
        // and additionnally reduce the value for convergence purposes (decreasecoeff) 
        	lame1 = dparam[8]*lame1*decreasecoeff;
   	        lame2 = dparam[8]*lame2*decreasecoeff;
            bVisc = dparam[8]*bVisc;
            activ = dparam[8] * deltat;
            facetag = 2;
            break;
    }

        // Determine which boundary does the lateral face belong to?
    int    boundarytag = -1; // for apical and basal should remain unchanged
    switch(key[1])
    {
        case 'T':
            boundarytag = 0;
            break;
        case 'B':
            boundarytag = 1;
            break;
        case 'L':
            boundarytag = 2;
            break;
        case 'R':
            boundarytag = 3;
            break;
    }

    //------------------------------------------------------------------
    // [2] OUTPUT definition
    std::vector<double>& rhs_k    = fillStr.Bk(0);   // <<=== MAIN OUTPUT
    std::vector<double>& g_rhs_k  = fillStr.Bk(1);   // <<=== MAIN OUTPUT
    std::vector<double>& Kmat_k   = fillStr.Ak(0,0); // <<=== MAIN OUTPUT
    std::vector<double>& Qmat_k   = fillStr.Ak(0,1); // <<=== MAIN OUTPUT
    std::vector<double>& Rmat_k   = fillStr.Ak(1,0); // <<=== MAIN OUTPUT
    std::vector<double>& Lmat_k   = fillStr.Ak(1,1); // <<=== MAIN OUTPUT

    //----------
    //Deformed
    double  x[3] = {};
    double e1[3] = {};
    double e2[3] = {};
    double  f[3] = {};

    for ( int i = 0; i < eNN; i++ )
    {
        for ( int n = 0; n < nDim; n++ )
        {
            x[n]  +=         p_k[i] * ue_nodes[DOF*i+n];
            e1[n] += dp_k[pDim*i+0] * ue_nodes[DOF*i+n];
            e2[n] += dp_k[pDim*i+1] * ue_nodes[DOF*i+n];

            f[n]  +=         p_k[i] * fe_nodes[DOF*i+n];
        }
    }


    //Metric tensor
    double Ip[4];
    Ip[0] = Math::Dot3D(e1,e1);
    Ip[1] = Math::Dot3D(e1,e2);
    Ip[2] = Ip[1];
    Ip[3] = Math::Dot3D(e2,e2);

    double det = Math::DetMat2x2(Ip);
    double J   = sqrt(det);

    //Inverse
    double iIp[4];
    Math::Invert2x2(iIp,Ip);

    double normal[3];
    Math::Cross(normal,e1,e2);
    Math::AX(normal,nDim,1.0/J);

    double xn = Math::Dot3D(x,normal);

    //----------
    //Previous time-step
    double  x0[3] = {};
    double e10[3] = {};
    double e20[3] = {};

    for ( int i = 0; i < eNN; i++ )
    {
        for ( int n = 0; n < nDim; n++ )
        {
            x0[n]  +=         p_k[i] * ue_nodes0[DOF*i+n];
            e10[n] += dp_k[pDim*i+0] * ue_nodes0[DOF*i+n];
            e20[n] += dp_k[pDim*i+1] * ue_nodes0[DOF*i+n];
        }
    }

    //Metric tensor
    double Ip0[4];
    Ip0[0] = Math::Dot3D(e10,e10);
    Ip0[1] = Math::Dot3D(e10,e20);
    Ip0[2] = Ip0[1];
    Ip0[3] = Math::Dot3D(e20,e20);

    double det0 = Math::DetMat2x2(Ip0);
    double J0   = sqrt(det0);

    //Inverse
    double iIp0[4];
    Math::Invert2x2(iIp0,Ip0);

    double normal0[3];
    Math::Cross(normal0,e10,e20);
    Math::AX(normal0,nDim,1.0/J0);

    double x0n0 = Math::Dot3D(x0,normal0);
    

    //----------
    //Reference
    double  xR[3] = {};
    double e1R[3] = {};
    double e2R[3] = {};

    for ( int i = 0; i < eNN; i++ )
    {
        for ( int n = 0; n < nDim; n++ )
        {
            xR[n]  +=         p_k[i] * xe_nodes[nDim*i+n];
            e1R[n] += dp_k[pDim*i+0] * xe_nodes[nDim*i+n];
            e2R[n] += dp_k[pDim*i+1] * xe_nodes[nDim*i+n];
        }
    }

    //Metric tensor
    double IpR[4];
    IpR[0] = Math::Dot3D(e1R,e1R);
    IpR[1] = Math::Dot3D(e1R,e2R);
    IpR[2] = IpR[1];
    IpR[3] = Math::Dot3D(e2R,e2R);

    double detR = Math::DetMat2x2(IpR);
    double JR   = sqrt(detR);

    //Inverse
    double iIpR[4];
    Math::Invert2x2(iIpR,IpR);

    double normalR[3];
    Math::Cross(normalR,e1R,e2R);
    Math::AX(normalR,nDim,1.0/JR);

    double xRnR = Math::Dot3D(xR,normalR);

    //----------
    //Rate-of-deformation
    double d[4]    = {};
    double d_Cc[4] = {};
    double d_CC[4] = {};

    Array::Copy(d,4,Ip);
    Math::AXPY(d,4, -1.0, Ip0);

    Math::MatProduct(d_Cc,2,2,2,iIp0,d);
    Math::MatProduct(d_CC,2,2,2,d_Cc,iIp0);

    double trd = (J-J0)/J0;

    double  v[3] = {};
    Array::Copy(v,3,x);
    Math::AXPY(v,3,-1.0,x0);

    double T[4]={};
    double iT[4]={};
    double GRef[4]={};
    double iGRef[4]={};
    double GRef0[4]={};
    double iGRef0[4]={};
    double dGRef[4]={};
    double Voigt[12]={}; //dG^ab/dG^cd
    double dG_GRef[12]={}, dG_GRef_Cc[12]={}, dG_GRef_cc[12]={};
    Voigt[4*0+2*0+0] = 1.0; Voigt[4*1+2*0+0] = 0.0; Voigt[4*2+2*0+0] = 0.0;
    Voigt[4*0+2*0+1] = 0.0; Voigt[4*1+2*0+1] = 1.0; Voigt[4*2+2*0+1] = 0.0;
    Voigt[4*0+2*1+0] = 0.0; Voigt[4*1+2*1+0] = 1.0; Voigt[4*2+2*1+0] = 0.0;
    Voigt[4*0+2*1+1] = 0.0; Voigt[4*1+2*1+1] = 0.0; Voigt[4*2+2*1+1] = 1.0;

    // basis vectors for each face
    double *g1e = &g_ae_nodes[0];
    double *g2e = &g_ae_nodes[3];

    // Change of basis from face to element
    LocMongeParam::matChangeBasis(T, g1e, g2e, e1R, e2R);
    Math::Invert2x2(iT,T);

    //Expressing G in element from that of the face (g_ue_nodes correspond to the value on the face) 
    double one = 1.0; // FE basis vector value
    CalcdGIGRef( dG_GRef, iT, one, Voigt );

    Math::AXPY(GRef,4,g_ue_nodes[0],&dG_GRef[4*0]);
    Math::AXPY(GRef,4,g_ue_nodes[1],&dG_GRef[4*1]);
    Math::AXPY(GRef,4,g_ue_nodes[2],&dG_GRef[4*2]);

    Math::AXPY(GRef0,4,g_ue_nodes0[0],&dG_GRef[4*0]);
    Math::AXPY(GRef0,4,g_ue_nodes0[1],&dG_GRef[4*1]);
    Math::AXPY(GRef0,4,g_ue_nodes0[2],&dG_GRef[4*2]);

    Math::Invert2x2(iGRef,GRef);
    Math::Invert2x2(iGRef0,GRef0);

    Array::Copy(dGRef,4,GRef);
    Math::AXPY(dGRef,4,-1.0,GRef0);

    Math::MatProduct(&dG_GRef_Cc[4*0],2,2,2,iGRef0,&dG_GRef[4*0]); // GRef0
    Math::MatProduct(&dG_GRef_Cc[4*1],2,2,2,iGRef0,&dG_GRef[4*1]);
    Math::MatProduct(&dG_GRef_Cc[4*2],2,2,2,iGRef0,&dG_GRef[4*2]);

    Math::MatProduct(&dG_GRef_cc[4*0],2,2,2,&dG_GRef_Cc[4*0],iGRef0);
    Math::MatProduct(&dG_GRef_cc[4*1],2,2,2,&dG_GRef_Cc[4*1],iGRef0);
    Math::MatProduct(&dG_GRef_cc[4*2],2,2,2,&dG_GRef_Cc[4*2],iGRef0);

    //Elasticity invariants
    double inva1 = Math::Dot(4,GRef,Ip);
    double sqrtdetG = sqrt(Math::DetMat2x2(GRef));
    double inva2 = sqrtdetG * J;
    double inva22 = inva2*inva2;

    //NeoHookean
    double ener              = lame1 * ( inva1/inva2 - 2.0 ) + lame2 * (inva2-1.0) * (inva2-1.0) ;
    double dinva1_Ener       = lame1 /inva2;
    double dinva2_Ener       = -lame1 * inva1/inva22 + 2.0 * lame2 * (inva2-1.0);
    double ddinva1_Ener      = 0.0;
    double ddinva2_Ener      = 2.0 * lame1 * inva1/(inva2*inva22) + 2.0 * lame2;
    double dinva1dinva2_Ener = -lame1 /inva22;

    //Free energy at the previous timestep
    double inva10    = Math::Dot(4,GRef0,Ip0);
    double sqrtdetG0 = sqrt(Math::DetMat2x2(GRef0));
    double inva20    = sqrtdetG0 * J0;
    double ener0     = lame1 * ( inva10/inva20 - 2.0 ) + lame2 * (inva20-1.0) * (inva20-1.0) ;

    // Values needed for contact
    double distsquared = fillStr.paramStr->x_aux[faceID]*fillStr.paramStr->x_aux[faceID] + fillStr.paramStr->y_aux[faceID]*fillStr.paramStr->y_aux[faceID] + fillStr.paramStr->z_aux[faceID]*fillStr.paramStr->z_aux[faceID];
    double distBetweenCellBary = sqrt(distsquared);

    int nPointsinCell = fillStr.paramStr->i_aux[faceID];
    double K = contact_force_coefficient * deltat;
    int distPow = -3;
    double offset = 0.8;
    double CoordDiffToOtherCellBary[] = {fillStr.paramStr->x_aux[faceID], fillStr.paramStr->y_aux[faceID], fillStr.paramStr->z_aux[faceID]};
    double domeScalingfactor = 1.0 / 3.0;

   
    /////////////////////////////////// FILLING THE RHS /////////////////////////////////// 

    for ( int i = 0; i < eNN; i++ ) {
     
        double dxI_Ip[4 * 3] = {}, dxI_Ip_Cc[4 * 3] = {}, dxI_Ip_CC[4 * 3] = {};
        MetricGradient(dxI_Ip, &dp_k[2 * i], e1, e2);

        Math::MatProduct(&dxI_Ip_Cc[0], 2, 2, 2, iIp0, &dxI_Ip[0]);
        Math::MatProduct(&dxI_Ip_Cc[4], 2, 2, 2, iIp0, &dxI_Ip[4]);
        Math::MatProduct(&dxI_Ip_Cc[8], 2, 2, 2, iIp0, &dxI_Ip[8]);

        Math::MatProduct(&dxI_Ip_CC[0], 2, 2, 2, &dxI_Ip_Cc[0], iIp0);
        Math::MatProduct(&dxI_Ip_CC[4], 2, 2, 2, &dxI_Ip_Cc[4], iIp0);
        Math::MatProduct(&dxI_Ip_CC[8], 2, 2, 2, &dxI_Ip_Cc[8], iIp0);

        double dxI_J[3] = {};
        JacobianGradient(dxI_J, &dp_k[2 * i], e1, e2, normal);

        double dxI_normal[3 * 3] = {};
        double xdxI_normal[3] = {};
        NormalGradient(dxI_normal, &dp_k[2 * i], e1, e2, normal, J, dxI_J);
        NormalGradientU(xdxI_normal, &dp_k[2 * i], e1, e2, normal, 1.0 / J, x);

        ///[1] Tension power
        Math::AXPY(&rhs_k[DOF * i], nDim, activ * dens0, dxI_J);

        ///[2] Cell volume constraint
        Math::AXPY(&rhs_k[DOF * i], nDim, deltat * pressure * J * p_k[i], normal);
        Math::AXPY(&rhs_k[DOF * i], nDim, deltat * pressure * J, xdxI_normal);
        Math::AXPY(&rhs_k[DOF * i], nDim, deltat * pressure * xn, dxI_J);

        ///[3] Elasticity
        // dxI_psi
        rhs_k[DOF * i + 0] += J0 * dens0 * dinva1_Ener * Math::Dot(4, GRef, &dxI_Ip[0]);
        rhs_k[DOF * i + 1] += J0 * dens0 * dinva1_Ener * Math::Dot(4, GRef, &dxI_Ip[4]);
        rhs_k[DOF * i + 2] += J0 * dens0 * dinva1_Ener * Math::Dot(4, GRef, &dxI_Ip[8]);

        rhs_k[DOF * i + 0] += J0 * dens0 * dinva2_Ener * sqrtdetG * dxI_J[0];
        rhs_k[DOF * i + 1] += J0 * dens0 * dinva2_Ener * sqrtdetG * dxI_J[1];
        rhs_k[DOF * i + 2] += J0 * dens0 * dinva2_Ener * sqrtdetG * dxI_J[2];

        ///[4] Surface friction dissipation
        // dxI_v
        rhs_k[DOF * i + 0] += J0 * eta_f * (x[0] - x0[0]) * p_k[i];
        rhs_k[DOF * i + 1] += J0 * eta_f * (x[1] - x0[1]) * p_k[i];
        rhs_k[DOF * i + 2] += J0 * eta_f * (x[2] - x0[2]) * p_k[i];

        
        /// [5] Luminal pressure, only on basal cell faces, working against cell pressure on basal surface
        // Implemented as a Lagrangian multiplier
        // Note the - sign as basal surface normal points towards the luminal cavity
        if (facetag == 1)
        {
            Math::AXPY(&rhs_k[DOF * i], nDim, -deltat * lumpressure * J * p_k[i] * domeScalingfactor, normal);
            Math::AXPY(&rhs_k[DOF * i], nDim, -deltat * lumpressure * J * domeScalingfactor, xdxI_normal);
            Math::AXPY(&rhs_k[DOF * i], nDim, -deltat * lumpressure * xn * domeScalingfactor, dxI_J);
        }
		
		 // [6] contact : potential = K/d
	    if (fillStr.paramStr->bparam[faceID] && subFill.kPt == 0)// && subFill.elemID == 0 && subFill.kPt== 0) //subFill.kPt== 0
        {

            rhs_k[DOF * i + 0] += K * CoordDiffToOtherCellBary[0] * pow(distBetweenCellBary - offset,distPow + 1) / distBetweenCellBary * p_k[i] / nPointsinCell  / subFill.wk_sample;
            rhs_k[DOF * i + 1] += K * CoordDiffToOtherCellBary[1] * pow(distBetweenCellBary - offset,distPow + 1) / distBetweenCellBary * p_k[i] / nPointsinCell  / subFill.wk_sample;
            rhs_k[DOF * i + 2] += K * CoordDiffToOtherCellBary[2] * pow(distBetweenCellBary - offset,distPow + 1) / distBetweenCellBary * p_k[i] / nPointsinCell  / subFill.wk_sample;

        }

	    
        /////////////////////////////////// FILLING K_MAT ///////////////////////////////////


        for (int j = 0; j < eNN; j++)
        {
            double Kmat_IJ[9] = {};

            double dxJ_Ip [4*3]={};
            MetricGradient(dxJ_Ip, &dp_k[2 * j], e1, e2);

            double dxJ_J[3]={};
            JacobianGradient(dxJ_J,&dp_k[2 * j], e1, e2, normal);

            double dxIdxJ_Ip[4]={};
            MetricHessian(dxIdxJ_Ip,&dp_k[2 * i], &dp_k[2 * j]);

            double dxJ_normal[9]={};
            double dxIdxJ_J[9]={};
            NormalGradient(dxJ_normal, &dp_k[2 * j], e1, e2, normal, J, dxJ_J);
            JacobianHessian(dxIdxJ_J, &dp_k[2 * i], &dp_k[2 * j], e1, e2, normal, dxJ_normal);

            double dxIdxJ_normal[27]={};
            NormalHessian(dxIdxJ_normal, &dp_k[2 * i], &dp_k[2 * j], J, dxI_J, dxJ_J, dxIdxJ_J, normal, dxI_normal, dxJ_normal);

            ///[1] Tension power
            for (int n = 0; n < nDim; n++)
                for (int m = 0; m < nDim; m++)
                    {
			            Kmat_IJ[DOF * n + m] += activ * dens0 * dxIdxJ_J[nDim * n + m];
		            }

            ///[2] Cell volume constraint
            double xdxJ_normal[3]={};
            NormalGradientU(xdxJ_normal,&dp_k[2 * j], e1, e2, normal, 1.0 / J, x);

            double xdxIdxJ_normal[9] = {};
            Math::AXPY(xdxIdxJ_normal, 9, x[0], &dxIdxJ_normal[0]);
            Math::AXPY(xdxIdxJ_normal, 9, x[1], &dxIdxJ_normal[9]);
            Math::AXPY(xdxIdxJ_normal, 9, x[2], &dxIdxJ_normal[18]);

            for (int n = 0; n < nDim; n++)
            {
                for (int m = 0; m < nDim; m++)
                {
                    Kmat_IJ[DOF * n + m] += deltat * pressure * dxJ_J[m] * p_k[i] * normal[n];
                    Kmat_IJ[DOF * n + m] += deltat * pressure * J * p_k[i] * dxJ_normal[n * nDim + m];

                    Kmat_IJ[DOF * n + m] += deltat * pressure * J * p_k[j] * dxI_normal[m * nDim + n];
                    Kmat_IJ[DOF * n + m] += deltat * pressure * xdxI_normal[n] * dxJ_J[m];
                    Kmat_IJ[DOF * n + m] += deltat * pressure * J * xdxIdxJ_normal[n * nDim + m];

                    Kmat_IJ[DOF * n + m] += deltat * pressure * xdxJ_normal[m] * dxI_J[n];
                    Kmat_IJ[DOF * n + m] += deltat * pressure * dxI_J[n] * p_k[j] * normal[m];
                    Kmat_IJ[DOF * n + m] += deltat * pressure * xn * dxIdxJ_J[n * nDim + m];
                }
            }

            ///[3] Elasticity and remodelling
            for (int n = 0; n < nDim; n++)
            {
                for (int m = 0; m < nDim; m++)
                {
                    // dxIdxJ_psi
                    Kmat_IJ[DOF * n + m] += J0 * dens0 * ddinva1_Ener * Math::Dot(4,GRef,&dxI_Ip[4*n])  * Math::Dot(4,GRef,&dxJ_Ip[4*m]);
                    Kmat_IJ[DOF * n + m] += J0 * dens0 * dinva1dinva2_Ener * Math::Dot(4,GRef,&dxI_Ip[4*n]) * sqrtdetG * dxJ_J[m];
                    Kmat_IJ[DOF * n + m] += J0 * dens0 * dinva2_Ener * sqrtdetG * dxIdxJ_J[n * nDim + m];
                    Kmat_IJ[DOF * n + m] += J0 * dens0 * ddinva2_Ener * sqrtdetG * dxI_J[n] * sqrtdetG * dxJ_J[m];
                    Kmat_IJ[DOF * n + m] += J0 * dens0 * dinva1dinva2_Ener * sqrtdetG * dxI_J[n] *  Math::Dot(4,GRef,&dxJ_Ip[4*m]);
                }

                // dxIdxJ_psi
                Kmat_IJ[DOF * n + n] += J0 * dens0 * dinva1_Ener * Math::Dot(4,GRef,dxIdxJ_Ip);
            }


            ///[4] Surface friction dissipation
            for (int n = 0; n < nDim; n++)
            // dxIdxJ_v
                Kmat_IJ[DOF * n + n] += J0 * eta_f * p_k[i] * p_k[j];


            /// [5] Luminal pressure
            if (facetag == 1)
            {
                for (int n = 0; n < nDim; n++)
                {
                    for (int m = 0; m < nDim; m++)
                    {
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * dxJ_J[m] * p_k[i] * normal[n] * domeScalingfactor;
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * J * p_k[i] * dxJ_normal[n * nDim + m] * domeScalingfactor;

                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * J * p_k[j] * dxI_normal[m * nDim + n] * domeScalingfactor;
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * xdxI_normal[n] * dxJ_J[m] * domeScalingfactor;
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * J * xdxIdxJ_normal[n * nDim + m] * domeScalingfactor;

                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * xdxJ_normal[m] * dxI_J[n] * domeScalingfactor;
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * dxI_J[n] * p_k[j] * normal[m] * domeScalingfactor;
                        Kmat_IJ[DOF * n + m] -= deltat * lumpressure * xn * dxIdxJ_J[n * nDim + m] * domeScalingfactor;
                    }
                }
            }

        /// [6] contact : potential = K/d
        if (fillStr.paramStr->bparam[faceID] && subFill.kPt == 0)
        {
		    for (int n = 0; n < nDim; n++) {
                        for (int m = 0; m < nDim; m++) {

			Kmat_IJ[DOF * n + m] += K * CoordDiffToOtherCellBary[n] * pow(distBetweenCellBary - offset,distPow + 1) * p_k[i] / nPointsinCell / 
                pow(distBetweenCellBary,distPow) / nPointsinCell * p_k[j] * CoordDiffToOtherCellBary[m]  / subFill.wk_sample;
            
            Kmat_IJ[DOF * n + m] += K / distBetweenCellBary / nPointsinCell * CoordDiffToOtherCellBary[n] * p_k[i] * 
                2 /  pow(distBetweenCellBary - offset,distPow) / distBetweenCellBary * p_k[j] /  nPointsinCell * CoordDiffToOtherCellBary[m] / subFill.wk_sample;

							}

            Kmat_IJ[DOF * n + n] -= K * pow(distBetweenCellBary - offset,distPow + 1) / distBetweenCellBary * p_k[i] / nPointsinCell *
                p_k[j] / nPointsinCell / subFill.wk_sample;

					}

        }
 
            // FILL term of K_e[i,j]
            for (int n = 0; n < DOF; n++)
                for (int m = 0; m < DOF; m++)
                        Kmat_k[(i * DOF + n) * N_k + j * DOF + m] += Kmat_IJ[n * DOF + m];


        }

    /////////////////////////////////// FILLING Q_MAT and R_MAT ///////////////////////////////////

        ///[2] Cell volume constraint
        Qmat_k[(DOF * i + 0) * gDOF + 3] += deltat * J * p_k[i] * normal[0];
        Qmat_k[(DOF * i + 1) * gDOF + 3] += deltat * J * p_k[i] * normal[1];
        Qmat_k[(DOF * i + 2) * gDOF + 3] += deltat * J * p_k[i] * normal[2];

        Qmat_k[(DOF * i + 0) * gDOF + 3] += deltat * J * xdxI_normal[0];
        Qmat_k[(DOF * i + 1) * gDOF + 3] += deltat * J * xdxI_normal[1];
        Qmat_k[(DOF * i + 2) * gDOF + 3] += deltat * J * xdxI_normal[2];

        Qmat_k[(DOF * i + 0) * gDOF + 3] += deltat * xn * dxI_J[0];
        Qmat_k[(DOF * i + 1) * gDOF + 3] += deltat * xn * dxI_J[1];
        Qmat_k[(DOF * i + 2) * gDOF + 3] += deltat * xn * dxI_J[2];


        Rmat_k[3 * N_k + DOF*i + 0] += deltat * J * p_k[i] * normal[0];
        Rmat_k[3 * N_k + DOF*i + 1] += deltat * J * p_k[i] * normal[1];
        Rmat_k[3 * N_k + DOF*i + 2] += deltat * J * p_k[i] * normal[2];

        Rmat_k[3 * N_k + DOF*i + 0] += deltat * J * xdxI_normal[0];
        Rmat_k[3 * N_k + DOF*i + 1] += deltat * J * xdxI_normal[1];
        Rmat_k[3 * N_k + DOF*i + 2] += deltat * J * xdxI_normal[2];

        Rmat_k[3 * N_k + DOF*i + 0] += deltat * xn * dxI_J[0];
        Rmat_k[3 * N_k + DOF*i + 1] += deltat * xn * dxI_J[1];
        Rmat_k[3 * N_k + DOF*i + 2] += deltat * xn * dxI_J[2];


        ///[3] Elasticity and remodelling
        for (int n = 0; n < nDim; n++)
        {
            for (int m = 0; m < nDim; m++)
            {
                // dxIdG_psi
                Qmat_k[(DOF * i + n) * gDOF + m] += J0 * dens0 * dinva1_Ener * Math::Dot(4,&dxI_Ip[4*n],&dG_GRef[4*m]);
                Qmat_k[(DOF * i + n) * gDOF + m] += J0 * dens0 * ddinva1_Ener * Math::Dot(4,GRef,&dxI_Ip[4*n]) * Math::Dot(4,Ip,&dG_GRef[4*m]);
                Qmat_k[(DOF * i + n) * gDOF + m] += J0 * dens0 * dinva1dinva2_Ener * Math::Dot(4,GRef,&dxI_Ip[4*n]) * inva2 * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*m]);
                Qmat_k[(DOF * i + n) * gDOF + m] += J0 * dens0 * (ddinva2_Ener * inva2 + dinva2_Ener) * sqrtdetG * dxI_J[n] * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*m]);
                Qmat_k[(DOF * i + n) * gDOF + m] += J0 * dens0 * dinva1dinva2_Ener * sqrtdetG * dxI_J[n]  * Math::Dot(4,Ip,&dG_GRef[4*m]);

                Rmat_k[n * N_k + DOF*i + m]    += J0 * dens0 * dinva1_Ener * Math::Dot(4,&dxI_Ip[4*m],&dG_GRef[4*n]);
                Rmat_k[n * N_k + DOF*i + m]    += J0 * dens0 * ddinva1_Ener * Math::Dot(4,GRef,&dxI_Ip[4*m]) * Math::Dot(4,Ip,&dG_GRef[4*n]);
                Rmat_k[n * N_k + DOF*i + m]    += J0 * dens0 * dinva1dinva2_Ener * Math::Dot(4,GRef,&dxI_Ip[4*m]) * inva2 * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*n]);
                Rmat_k[n * N_k + DOF*i + m]    += J0 * dens0 * (ddinva2_Ener * inva2 + dinva2_Ener) * sqrtdetG * dxI_J[m] * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*n]);
                Rmat_k[n * N_k + DOF*i + m]    += J0 * dens0 * dinva1dinva2_Ener * sqrtdetG * dxI_J[m]  * Math::Dot(4,Ip,&dG_GRef[4*n]);

            }
        }

        //[7] Density
        for (int n = 0; n < nDim; n++)
            Rmat_k[4 * N_k + DOF*i + n] += dens * dxI_J[n];

        // [5] Luminal pressure
        if (facetag == 1)
        {
            Qmat_k[(DOF * i + 0) * gDOF + 5] -= deltat * J * p_k[i] * normal[0] * domeScalingfactor;
            Qmat_k[(DOF * i + 1) * gDOF + 5] -= deltat * J * p_k[i] * normal[1] * domeScalingfactor;
            Qmat_k[(DOF * i + 2) * gDOF + 5] -= deltat * J * p_k[i] * normal[2] * domeScalingfactor;

            Qmat_k[(DOF * i + 0) * gDOF + 5] -= deltat * J * xdxI_normal[0] * domeScalingfactor;
            Qmat_k[(DOF * i + 1) * gDOF + 5] -= deltat * J * xdxI_normal[1] * domeScalingfactor;
            Qmat_k[(DOF * i + 2) * gDOF + 5] -= deltat * J * xdxI_normal[2] * domeScalingfactor;

            Qmat_k[(DOF * i + 0) * gDOF + 5] -= deltat * xn * dxI_J[0] * domeScalingfactor;
            Qmat_k[(DOF * i + 1) * gDOF + 5] -= deltat * xn * dxI_J[1] * domeScalingfactor;
            Qmat_k[(DOF * i + 2) * gDOF + 5] -= deltat * xn * dxI_J[2] * domeScalingfactor;

            Rmat_k[5 * N_k + DOF*i + 0] -= deltat * J * p_k[i] * normal[0] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 1] -= deltat * J * p_k[i] * normal[1] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 2] -= deltat * J * p_k[i] * normal[2] * domeScalingfactor;

            Rmat_k[5 * N_k + DOF*i + 0] -= deltat * J * xdxI_normal[0] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 1] -= deltat * J * xdxI_normal[1] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 2] -= deltat * J * xdxI_normal[2] * domeScalingfactor;

            Rmat_k[5 * N_k + DOF*i + 0] -= deltat * xn * dxI_J[0] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 1] -= deltat * xn * dxI_J[1] * domeScalingfactor;
            Rmat_k[5 * N_k + DOF*i + 2] -= deltat * xn * dxI_J[2] * domeScalingfactor;
        }
    }


    /////////////////////////////////// FILLING G_RHS and L_MAT/////////////////////////////////// 

    ///[2] Cell Volume constraint
    g_rhs_k[3] += deltat * ( xn*J - x0n0*J0);

    // [5] Luminal pressure, deltat * deltat * lumpressure * (Vdome - Vdomeincreaserate)
    if (facetag == 1)
        g_rhs_k[5] -= deltat * ( xn*J - x0n0*J0 ) * domeScalingfactor;

    // constant term corresponding to rate of luminal volume increase. Added only once
     if (faceID == 0 && subFill.elemID==0 && subFill.kPt== 0)
       g_rhs_k[5] -= deltat *  deltat * vdome / subFill.wk_sample;

    ///[7] Density
    g_rhs_k[4] += J0 * ( dens*J/J0 - dens0 - kp * conc + kd * dens );

    Lmat_k[gDOF * 4 + 4] += J + kd * J0;

    double dG_iGRef[12], MauxI[4];
    Math::MatProduct(MauxI,2,2,2,iGRef,&dG_GRef[4*0]);
    Math::MatProduct(&dG_iGRef[4*0],2,2,2,MauxI,iGRef);
    Math::MatProduct(MauxI,2,2,2,iGRef,&dG_GRef[4*1]);
    Math::MatProduct(&dG_iGRef[4*1],2,2,2,MauxI,iGRef);
    Math::MatProduct(MauxI,2,2,2,iGRef,&dG_GRef[4*2]);
    Math::MatProduct(&dG_iGRef[4*2],2,2,2,MauxI,iGRef);
    Math::AX(dG_iGRef,12,-1.0);

    /// [3] Elasticity
    // dG_psi
    g_rhs_k[0] += J0 * dens0 * dinva1_Ener * Math::Dot(4,Ip,&dG_GRef[4*0]);
    g_rhs_k[1] += J0 * dens0 * dinva1_Ener * Math::Dot(4,Ip,&dG_GRef[4*1]);
    g_rhs_k[2] += J0 * dens0 * dinva1_Ener * Math::Dot(4,Ip,&dG_GRef[4*2]);

    g_rhs_k[0] += J0 * dens0 * dinva2_Ener * inva2 * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*0]);
    g_rhs_k[1] += J0 * dens0 * dinva2_Ener * inva2 * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*1]);
    g_rhs_k[2] += J0 * dens0 * dinva2_Ener * inva2 * 0.5 * Math::Dot(4,iGRef,&dG_GRef[4*2]);

    for (int n = 0; n < nDim; n++)
    {
        for (int m = 0; m < nDim; m++)
        {
            Lmat_k[gDOF * n + m ] += J0 * dens0 * ddinva1_Ener * Math::Dot(4, Ip, &dG_GRef[4 * n]) * Math::Dot(4, Ip, &dG_GRef[4 * m]);
            Lmat_k[gDOF * n  + m ] += J0 * dens0 * dinva1dinva2_Ener * Math::Dot(4, Ip, &dG_GRef[4 * n]) * inva2 * 0.5 *  Math::Dot(4, iGRef, &dG_GRef[4 * m]);
            Lmat_k[gDOF * n  + m ] += J0 * dens0 * dinva2_Ener * inva2 * 0.5 * Math::Dot(4, &dG_GRef[4 * n], &dG_iGRef[4 * m]);
            Lmat_k[gDOF * n  + m ] += J0 * dens0 * (ddinva2_Ener * inva2 + dinva2_Ener) * 0.5 * Math::Dot(4, iGRef, &dG_GRef[4 * n]) * inva2 * 0.5 * Math::Dot(4, iGRef, &dG_GRef[4 * m]);
            Lmat_k[gDOF * n  + m ] += J0 * dens0 * dinva1dinva2_Ener * inva2 * 0.5 * Math::Dot(4, iGRef, &dG_GRef[4 * n]) * Math::Dot(4, Ip, &dG_GRef[4 * m]);
        }

    }

    ///[3] Remodeling
    // dGId_Dremod
    g_rhs_k[0] += bVisc * J0 * dens0 * Math::Dot(4,dGRef,&dG_GRef_cc[4*0]);
    g_rhs_k[1] += bVisc * J0 * dens0 * Math::Dot(4,dGRef,&dG_GRef_cc[4*1]);
    g_rhs_k[2] += bVisc * J0 * dens0 * Math::Dot(4,dGRef,&dG_GRef_cc[4*2]);

    // dGIdGJdxJ_Dremod
    for (int n = 0; n < nDim; n++)
        for (int m = 0; m < nDim; m++)
			Lmat_k[gDOF * n + m] += bVisc * J0 * dens0 * Math::Dot(4,&dG_GRef_cc[4*n],&dG_GRef[4*m]);



}


double Getloadstep(double t, double inflation_begin, double inflation_end,
                   double deflation_begin, double deflation_end, double deflation_mag, double inflation_mag)
{
    double incr = 0.0;

    // pull then push
    if (t > inflation_begin and t < inflation_end)
    {
        incr = (t - inflation_begin) * inflation_mag/(inflation_end - inflation_begin);
    }
    else if (t >= inflation_end and t <= deflation_begin)
    {
        incr = inflation_mag;
    }
    else if (t > deflation_begin and t < deflation_end)
    {
        incr = inflation_mag + (t - deflation_begin) * -(deflation_mag + inflation_mag)/(deflation_end - deflation_begin);
    }
    else if (t >= deflation_end)
     {
        incr = -deflation_mag;
    }
    else
    {
        incr = 0.0;
    }

    return incr;
}

