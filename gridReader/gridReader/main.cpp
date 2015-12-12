//
//  main.cpp
//  gridReader
//
//  Created by Franz Neubert on 06/10/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <math.h>

#include <iostream>
#include <vtkSmartPointer.h>

#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridGeometryFilter.h>

#include <vtkHexahedron.h>
#include <vtkTetra.h>

#include <vtkEdgeTable.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>

#include <vtkIdList.h>
#include <vtkIdTypeArray.h>

#include <vtkAlgorithmOutput.h>

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>

#include <vtkDelaunay3D.h>

#include <vtkGeometryFilter.h>
#include <vtkDataSetSurfaceFilter.h>

#include <vtkAppendFilter.h>

#include <vtkDataSetTriangleFilter.h>

#include <vtkExtractEdges.h>

#include <vtkFeatureEdges.h>

#include <vtkFieldDataToAttributeDataFilter.h>

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkInformation.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "Calculator.h"
#include "CostCalculator.hpp"
#include "EdgeCollapse.hpp"
#include "Helper.hpp"
#include "Vertex.hpp"

struct CompareCost {
    bool operator()(EdgeCollapse &col1, EdgeCollapse &col2) {
        return col1.getCost() > col2.getCost();
    }
};



void hexahedronToTetrahedronGrid(vtkUnstructuredGrid *hexagrid, vtkUnstructuredGrid *tetraGrid){
    vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangleFilter->SetInputData(hexagrid);
    triangleFilter->Update();
    std::cout << "Hexagrid with  " << hexagrid->GetNumberOfCells() << "  cells to \n";
    std::cout << "Tetragrid with " << triangleFilter->GetOutput()->GetNumberOfCells() << " cells" << std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    tetraGrid = grid;
}


/*
 Saves a unstructured grid to a given .vtu file (full path required).
 */
void writeUgrid(vtkUnstructuredGrid *ugrid, const char *filename) {
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(ugrid);
    writer->Write();
}

void startRenderingAlgo(vtkAlgorithmOutput* grid) {
    grid->PrintSelf(std::cout, *vtkIndent::New());
    vtkSmartPointer<vtkExtractEdges> extract = vtkSmartPointer<vtkExtractEdges>::New();
    extract->SetInputConnection(grid);
    extract->Update();
    
    vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geoFilter->SetInputConnection(grid);
    geoFilter->Update();
    
    vtkSmartPointer<vtkFeatureEdges> features = vtkSmartPointer<vtkFeatureEdges>::New();
    features->SetInputConnection(geoFilter->GetOutputPort());
    features->BoundaryEdgesOn();
    features->FeatureEdgesOn();
    features->ManifoldEdgesOn();
    features->NonManifoldEdgesOn();
    features->Update();
    
    std::cout << features->GetOutput()->GetNumberOfCells() << std::endl;
    
    vtkSmartPointer<vtkPolyDataMapper> featureMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    featureMapper->SetInputConnection(features->GetOutputPort());
    
    vtkSmartPointer<vtkActor> featureActor = vtkSmartPointer<vtkActor>::New();
    featureActor->SetMapper(featureMapper);
    
    vtkSmartPointer<vtkPolyDataMapper> edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    edgeMapper->SetInputConnection(extract->GetOutputPort());
    
    vtkSmartPointer<vtkActor> edgeActor = vtkSmartPointer<vtkActor>::New();
    edgeActor->SetMapper(edgeMapper);
    edgeActor->GetProperty()->SetColor(0, 0, 0);
    
    //Mapper
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(grid);
    
    //Actor
    vtkSmartPointer<vtkActor> ugridActor = vtkSmartPointer<vtkActor>::New();
    ugridActor->SetMapper(mapper);
    ugridActor->GetProperty()->SetOpacity(0.5);
    
    
    //Renderer
    vtkSmartPointer<vtkRenderer> ugridRenderer = vtkSmartPointer<vtkRenderer>::New();
    ugridRenderer->AddActor(featureActor);
    ugridRenderer->AddActor(ugridActor);
    //ugridRenderer->AddActor(edgeActor);
    ugridRenderer->SetBackground(0.1, 0.2, 0.4);
    
    //Render Window & Interactor
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ugridRenderer);
    renWin->SetSize(1024, 1024);
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    iren->Initialize();
    iren->Start();
}

void startRenderingGrid(vtkUnstructuredGrid *grid) {
    vtkSmartPointer<vtkExtractEdges> extract = vtkSmartPointer<vtkExtractEdges>::New();
    extract->SetInputData(grid);
    extract->Update();
    
    vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geoFilter->SetInputData(grid);
    geoFilter->Update();
    
    vtkSmartPointer<vtkFeatureEdges> features = vtkSmartPointer<vtkFeatureEdges>::New();
    features->SetInputConnection(geoFilter->GetOutputPort());
    features->BoundaryEdgesOn();
    features->FeatureEdgesOn();
    features->ManifoldEdgesOn();
    features->NonManifoldEdgesOn();
    features->Update();
    
    vtkSmartPointer<vtkPolyDataMapper> featureMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    featureMapper->SetInputConnection(features->GetOutputPort());
    
    vtkSmartPointer<vtkActor> featureActor = vtkSmartPointer<vtkActor>::New();
    featureActor->SetMapper(featureMapper);
    
    vtkSmartPointer<vtkPolyDataMapper> edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    edgeMapper->SetInputConnection(extract->GetOutputPort());
    
    vtkSmartPointer<vtkActor> edgeActor = vtkSmartPointer<vtkActor>::New();
    edgeActor->SetMapper(edgeMapper);
    edgeActor->GetProperty()->SetColor(0, 0, 0);
    
    //Mapper
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(grid);
    
    //Actor
    vtkSmartPointer<vtkActor> ugridActor = vtkSmartPointer<vtkActor>::New();
    ugridActor->SetMapper(mapper);
    ugridActor->GetProperty()->SetOpacity(0.5);
    
    
    //Renderer
    vtkSmartPointer<vtkRenderer> ugridRenderer = vtkSmartPointer<vtkRenderer>::New();
    // ugridRenderer->AddActor(featureActor);
    ugridRenderer->AddActor(ugridActor);
    ugridRenderer->AddActor(edgeActor);
    ugridRenderer->SetBackground(0.1, 0.2, 0.4);
    
    //Render Window & Interactor
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ugridRenderer);
    renWin->SetSize(1500, 1024);
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    iren->Initialize();
    iren->Start();
}

vtkSmartPointer<vtkUnstructuredGrid> generateSimpleGrid() {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int k = 0; k < 10; k++) {
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                points->InsertNextPoint(i, j, k);
            }
        }
    }
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetPoints(points);
    grid->SetDimensions(10, 10, 10);
    
    vtkSmartPointer<vtkAppendFilter> appender = vtkSmartPointer<vtkAppendFilter>::New();
    appender->SetInputData(grid);
    appender->Update();
    
    return appender->GetOutput();
}

bool isInnerVertex(double solidAngle) {
    if (solidAngle >= 3.1) {
        return true;
    }
    return false;
}

bool isCorner(double solidAngle) {
    if (solidAngle <= M_PI/2) {
        return true;
    }
    return false;
}

bool isBoundaryEdge(vtkIdType pointA, vtkIdType pointB, std::set<vtkIdType> *surfaceIds) {
    if (surfaceIds->find(pointA) != surfaceIds->end() && surfaceIds->find(pointB) != surfaceIds->end()) {
        // std::cout << pointA << "-->" << pointB << " is BOUNDARY" << std::endl;
        return true;
    }
    return false;
}

bool isInteriorEdge(vtkIdType pointA, vtkIdType pointB, std::set<vtkIdType> *surfaceIds) {
    if (surfaceIds->find(pointA) == surfaceIds->end() && surfaceIds->find(pointB) == surfaceIds->end()) {
        // std::cout << pointA << "-->" << pointB << " is INTERIOR" << std::endl;
        return true;
    }
    return false;
}

std::set<vtkIdType> getPointNeighbours(vtkUnstructuredGrid *grid, vtkIdType seed) {
    vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
    grid->GetPointCells(seed, pointCells);
    
    std::set<vtkIdType> pointCellIds;
    for (vtkIdType point = 0; point < pointCells->GetNumberOfIds(); point++) {
        pointCellIds.insert(pointCells->GetId(point));
    }
    
    
    return pointCellIds;
}

double getVertexSolidAngle(vtkUnstructuredGrid *grid, std::set<vtkIdType> cells, vtkIdType seed) {
    
    vtkPoints *points = grid->GetPoints();
    double solidAngleSum = 0;
    for (auto cell : cells) {
        vtkIdList* pointIds = grid->GetCell(cell)->GetPointIds();
        std::list<vtkIdType> cellPointIds = Helper::toStdList(pointIds);
        double pointCoords[3*3];
        double seedCoords[3];
        points->GetPoint(seed, seedCoords);
        int i = 0;
        for (auto point : cellPointIds) {
            if (point != seed) {
                points->GetPoint(point, &pointCoords[i]);
                i += 3;
            }
        }
        double a[3];
        double b[3];
        double c[3];
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[0], a);
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[3], b);
        Calculator::calcVectorBetweenPoints(seedCoords, &pointCoords[6], c);
        solidAngleSum += Calculator::calcSolidAngle(a, b, c);
    }
    
    return solidAngleSum;
}

std::set<vtkIdType> getSurfaceVertexIds(vtkUnstructuredGrid *grid) {
    std::set<vtkIdType> vertexIds;
    
    vtkPoints *gridPoints = grid->GetPoints();
    for (vtkIdType seed = 0; seed < gridPoints->GetNumberOfPoints(); seed++) {
        double seedCoords[3];
        gridPoints->GetPoint(seed, seedCoords);
        std::set<vtkIdType> pointCells = getPointNeighbours(grid, seed);
        double angle = getVertexSolidAngle(grid, pointCells, seed);
        if (isCorner(angle)) {
            // std::cout << seed << " is corner" << std::endl;
            vertexIds.insert(seed);
        } else if (!isInnerVertex(angle)) {
            // std::cout << seed << " is boundary" << std::endl;
            vertexIds.insert(seed);
        } else {
            // std::cout << seed << " is inner" << std::endl;
        }
    }
    
    return vertexIds;
}

void recalcNcells(std::vector<EdgeCollapse> *prio_q, EdgeCollapse *collapse, vtkUnstructuredGrid *tetraGrid) {
    
}


/*!
 Does the cost calucation step and builds up the edge collapse heap in the process
 \param tetraGrid The unstructured grid to calculate the cost and edged collapses from
 \returns The priority queue containing all possible edge collapses
 */
std::vector<EdgeCollapse>recalculateCosts(vtkUnstructuredGrid *tetraGrid) {
    std::vector<EdgeCollapse> prio_q;
    
    // detect ids of boundary points
    std::set<vtkIdType> surfacePoints = getSurfaceVertexIds(tetraGrid);
    
    vtkSmartPointer<vtkEdgeTable> edges = vtkSmartPointer<vtkEdgeTable>::New();
    
    // Iterate over tetrahedras
    for(vtkIdType tetra = 0; tetra < tetraGrid->GetNumberOfCells(); tetra++){
        
        if ((int)tetra%1000 == 0) {
            std::cout << "c = " << tetra << std::endl;
        }
        
        // Iterate over cell edges
        for (vtkIdType e = 0; e < tetraGrid->GetCell(tetra)->GetNumberOfEdges(); e++) {
            
            vtkIdType pointA = tetraGrid->GetCell(tetra)->GetEdge((int) e)->GetPointId(0);
            vtkIdType pointB = tetraGrid->GetCell(tetra)->GetEdge((int) e)->GetPointId(1);
            
            //check if edge was already traversed
            if (edges->IsEdge(pointA, pointB) == -1) {
                // mark edge as traversed
                edges->InsertEdge(pointA, pointB);
                
                // check if edge is boundary or interior (don´t collapse edges which have an interior and a boundary vertex)
                if ((isBoundaryEdge(pointA, pointB, &surfacePoints) || isInteriorEdge(pointA, pointB, &surfacePoints))) {
                    std::set<vtkIdType> icells = EdgeCollapse::getIcells(pointA, pointB, tetraGrid);
                    std::set<vtkIdType> ncells = EdgeCollapse::getNCells(pointA, pointB, tetraGrid);
                    double volumeCost = CostCalculator::calcVolumeCost(pointA, pointB, 500, &ncells, &icells , tetraGrid);
                    // double scalarCost = CostCalculator::calcScalarCost(pointA, pointB, 1, tetraGrid);
                    double collapseCost = volumeCost;
                    EdgeCollapse collapse = *new EdgeCollapse(pointA, pointB, collapseCost);
                    //collapse.setIcells(icells);
                    //collapse.setNcells(ncells);
                    prio_q.push_back(collapse);
                }
            }
        }
    }
    
    std::make_heap(prio_q.begin(), prio_q.end(), CompareCost());
    return prio_q;
}

/*!
 Does one edge collapse. Changes the coords of ncells accordingly and removes icells from the grid.
 \param collapse The EdgeCollapse to be executed.
 \param collapsingGrid The unstructured grid the collapse should be executed on.
 \returns The unstructured grid after the collapse.
 */
vtkSmartPointer<vtkUnstructuredGrid> doCollapse(EdgeCollapse *collapse, vtkUnstructuredGrid *collapsingGrid) {
    
    vtkIdType p1 = collapse->getPointA();
    vtkIdType p2 = collapse->getPointB();
    collapse->setIcells(EdgeCollapse::getIcells(p1, p2, collapsingGrid));
    collapse->setNcells(EdgeCollapse::getNCells(p1, p2, collapsingGrid));
    
    std::cout << "collapse: " << p1 << "--" << p2 << std::endl;
    
    

//  ------- NCELLS --------
//  Changes the coords of p1 to the midpoint between p1 and p2.
//  Replaces p2(id) with p1(id) in every ncell (which uses p2).

    double midpoint[3];
    double coordsP1[3];
    double coordsP2[3];
    collapsingGrid->GetPoint(p1, coordsP1);
    collapsingGrid->GetPoint(p2, coordsP2);
    Calculator::calcMidPoint(coordsP1, coordsP2, midpoint);
    collapsingGrid->GetPoints()->SetPoint(p1, midpoint);
    for (auto ncell : collapse->ncells) {
        vtkIdList *ncellIds = collapsingGrid->GetCell(ncell)->GetPointIds();
        if (ncellIds->IsId(p2) != -1) {
            ncellIds->SetId(ncellIds->IsId(p2), p1);
            collapsingGrid->ReplaceCell(ncell, (int)ncellIds->GetNumberOfIds(), ncellIds->GetPointer(0));
        }
    }
//  ---- END NCELLS ----
    

//  ------ ICELLS ------
//  Builds a vtkSelection from icell ids.
//  Inverts Selection and extracts it, retrieving all cells except the icells.

    vtkSmartPointer<vtkIdTypeArray> icellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    vtkSmartPointer<vtkExtractSelection> extract = vtkSmartPointer<vtkExtractSelection>::New();
    icellIds->SetNumberOfComponents(1);
    for (auto icell : collapse->icells) {
        icellIds->InsertNextValue(icell);
        std::cout << icell << std::endl;
    }
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);
    selectionNode->SetSelectionList(icellIds);
    selection->AddNode(selectionNode);
    extract->SetInputData(0, collapsingGrid);
    extract->SetInputData(1, selection);
    extract->Update();
    vtkSmartPointer<vtkUnstructuredGrid> newGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    newGrid->ShallowCopy(extract->GetOutput());
    /* ----- END ICELLS ----- */
    return newGrid;
}

int main(int argc, const char * argv[]) {
    
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    /*std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
     reader->SetFileName(inputFilename.c_str());
     reader->Update();
    */
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    
    if (/*reader->IsFileUnstructuredGrid()*/ true) {
        
        std::cout << "input file is unstructured grid" << std::endl;
        
        //Initialize UnstructuredGridReader
//        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
//        gridReader->SetFileName(inputFilename.c_str());
//        gridReader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> simpleGrid = generateSimpleGrid();
        std::cout << "_original_ number of cells: " << simpleGrid->GetNumberOfCells() << std::endl;
        
        // -----------------
        // quad to tetra
        // -----------------
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        //triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->SetInputData(simpleGrid);
        triangleFilter->Update();
        std::cout << "_after triangulation_ number of cells: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();

        
        // --------------------
        // Initiate algorithm
        // --------------------
        std::vector<EdgeCollapse> prio_q = recalculateCosts(tetraGrid);
        std::cout << "prio queue: " << prio_q.size() << std::endl;
        
        // start iteration over edge collapses
        for (int i = 0; i < 10; i++) {
            
            // report progress
            if (i%100 == 0) {
                std::cout << "collapses done: " << i << std::endl;
            }
            
            if (!prio_q.empty()) {
                
                // get edge collapse with lowest costs
                EdgeCollapse top = prio_q.front();
                
                // execute collapse
                tetraGrid = doCollapse(&top, tetraGrid);
                std::cout << "_after first simplification_: " << tetraGrid->GetNumberOfCells() << std::endl;
                
                // remove executed collapse (and keep heap property)
                // prio_q.erase(prio_q.begin());
                // std::make_heap(prio_q.begin(), prio_q.end(), CompareCost());

            }
        }
      
        writeUgrid(tetraGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20151212_with-ncell.vtu");
        // startRenderingGrid(collapsingGrid);
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


