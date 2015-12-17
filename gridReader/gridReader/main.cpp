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
#include <thread>
#include <mutex>

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

std::set<vtkIdType> removedPointIds;


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

bool isBoundary(double solidAngle) {
    if (solidAngle < 4 * M_PI) {
        return true;
    }
    return false;
}

bool isCorner(double solidAngle) {
    if (solidAngle <= M_PI/2 ||  4 * M_PI - solidAngle <= M_PI/2) {
        return true;
    }
    return false;
}

bool isCurveCorner(double solidAngle) {
    if ((M_PI/2 < solidAngle && solidAngle <= (3 * M_PI) / 2) ||
        ((M_PI/2 < 4 * M_PI - solidAngle) && (4 * M_PI - solidAngle <= (3 * M_PI) / 2))) {
        return true;
    }
    return false;
}

bool isBoundaryEdge(vtkIdType pointA, vtkIdType pointB, std::set<vtkIdType> *surfaceIds) {
    if (surfaceIds->find(pointA) != surfaceIds->end() && surfaceIds->find(pointB) != surfaceIds->end()) {
        return true;
    }
    return false;
}

bool isInteriorEdge(vtkIdType pointA, vtkIdType pointB, std::set<vtkIdType> *surfaceIds) {
    if (surfaceIds->find(pointA) == surfaceIds->end() && surfaceIds->find(pointB) == surfaceIds->end()) {
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
    pointCells = NULL;
    return pointCellIds;
}

/*!
 Calculates the sum of the solid angles of tetrahedron incident on a vertex (seed).
 \param grid The unstructured grid containing the tetrahedrons and the seed vertex
 \param seedCells A std::set of cell ids from cells incident on the seed vertex
 \param seed The vertex id of the seed vertex
 \return the sum of solid angles of all tetrahedrons incident on the seed (in radians/steridians e.g. pi/2)
 */
double getVertexSolidAngle(vtkUnstructuredGrid *grid, std::set<vtkIdType> seedCells, vtkIdType seed) {
    vtkPoints *points = grid->GetPoints();
    double solidAngleSum = 0;
    for (auto cell : seedCells) {
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
        solidAngleSum += Calculator::calcSolidAngle(a, b, c, seedCoords);
    }
    return solidAngleSum;
}

/*!
 Retrieves a list of vertices which sit on the boundary of the grid (not corner, outer-edge).
 \param grid The unstructured grid to retrieve the vertices from
 \return A std::set of vtkIdTypes, representing the ids of general boundary vertices
 */
std::set<vtkIdType> getSurfaceVertexIds(vtkUnstructuredGrid *grid) {
    std::set<vtkIdType> vertexIds;
    
    vtkPoints *gridPoints = grid->GetPoints();
    for (vtkIdType seed = 0; seed < gridPoints->GetNumberOfPoints(); seed++) {
        if (removedPointIds.find(seed) == removedPointIds.end()) {
            double seedCoords[3];
            gridPoints->GetPoint(seed, seedCoords);
            std::set<vtkIdType> pointCells = Helper::cellsUsingVertex(seed, grid);
            double angle = getVertexSolidAngle(grid, pointCells, seed);
            if (!isCorner(angle) && !isCurveCorner(angle) && isBoundary(angle)) {
                vertexIds.insert(seed);
            }
        }
    }
    return vertexIds;
}


void iterateTetras(vtkUnstructuredGrid *tetraGrid,
                   std::vector<EdgeCollapse> *prio_q,
                   std::set<vtkIdType> surfacePoints,
                   vtkIdType iterationStart,
                   vtkIdType iterationEnd) {
    
    std::set<vtkIdType> collapsePointIds;
    
    for(vtkIdType tetra = iterationStart; tetra < iterationEnd; tetra++){
        
        if ((int)tetra%1000 == 0) {
            std::cout << "c = " << tetra << std::endl;
        }
        
        // Iterate over cell edges
        for (vtkIdType e = 0; e < tetraGrid->GetCell(tetra)->GetNumberOfEdges(); e++) {
            
            vtkIdType pointA = tetraGrid->GetCell(tetra)->GetEdge((int) e)->GetPointId(0);
            vtkIdType pointB = tetraGrid->GetCell(tetra)->GetEdge((int) e)->GetPointId(1);
            
            //check if edge was already traversed
            if (collapsePointIds.find(pointA) == collapsePointIds.end() &&
                collapsePointIds.find(pointB) == collapsePointIds.end() &&
                removedPointIds.find(pointA) == removedPointIds.end() &&
                removedPointIds.find(pointB) == removedPointIds.end()) {
                // mark edge as traversed
                collapsePointIds.insert(pointA);
                collapsePointIds.insert(pointB);
                
                // check if edge is boundary or interior (don´t collapse edges which have an interior and a boundary vertex)
                if ((isBoundaryEdge(pointA, pointB, &surfacePoints) || isInteriorEdge(pointA, pointB, &surfacePoints))) {
                    std::set<vtkIdType> icells = EdgeCollapse::getIcells(pointA, pointB, tetraGrid);
                    std::set<vtkIdType> ncells = EdgeCollapse::getNCells(pointA, pointB, tetraGrid);
                    double volumeCost = CostCalculator::calcVolumeCost(pointA, pointB, 500, &ncells, &icells , tetraGrid);
                    // double scalarCost = CostCalculator::calcScalarCost(pointA, pointB, 250, tetraGrid);
                    double edgeLengthCost = CostCalculator::calcEdgeLengthCost(pointA, pointB, 5, tetraGrid);
                    double collapseCost = volumeCost + edgeLengthCost;
                    // EdgeCollapse *collapse = new EdgeCollapse(pointA, pointB, collapseCost);
                    prio_q->push_back(*new EdgeCollapse(pointA, pointB, collapseCost));
                }
            }
        }
    }
    
}

/*!
 Does the cost calucation step and builds up the edge collapse heap in the process
 \param tetraGrid The unstructured grid to calculate the cost and edged collapses from
 \returns The priority queue containing all possible edge collapses
 */
std::vector<EdgeCollapse>recalculateCosts(vtkUnstructuredGrid *tetraGrid) {
    std::vector<EdgeCollapse> prio_q;
    std::set<vtkIdType> collapsePointIds;
    
    // detect ids of boundary points
    std::set<vtkIdType> surfacePoints = getSurfaceVertexIds(tetraGrid);
    iterateTetras(tetraGrid, &prio_q, surfacePoints, 0, tetraGrid->GetNumberOfCells());
    
    std::make_heap(prio_q.begin(), prio_q.end(), CompareCost());
    std::cout << "prio queue: " << prio_q.size() << std::endl;
    return prio_q;
}

/*!
 Does one edge collapse. Changes the coords of ncells accordingly and removes icells from the grid.
 \param collapse The EdgeCollapse to be executed.
 \param collapsingGrid The unstructured grid the collapse should be executed on.
 \returns The unstructured grid after the collapse.
 */
vtkSmartPointer<vtkUnstructuredGrid> doCollapse(EdgeCollapse *collapse, vtkUnstructuredGrid *collapsingGrid) {
    vtkSmartPointer<vtkUnstructuredGrid> newGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkIdType p1 = collapse->getPointA();
    vtkIdType p2 = collapse->getPointB();
    collapse->setNcells(EdgeCollapse::getNCells(p1, p2, collapsingGrid));
    collapse->setIcells(EdgeCollapse::getIcells(p1, p2, collapsingGrid));
    
    std::cout << "collapse: " << p1 << "--" << p2 << std::endl;
    std::cout << collapsingGrid->GetNumberOfPoints() << std::endl;
    
    //  ------- NCELLS --------
    //  Changes the coords of p1 to the midpoint between p1 and p2.
    //  Replaces p2(id) with p1(id) in every ncell (using p2).
    //  Removes the unused vertex p2 via inverse ExtracSelection
    double midpoint[3];
    double coordsP1[3];
    double coordsP2[3];
    collapsingGrid->GetPoint(p1, coordsP1);
    collapsingGrid->GetPoint(p2, coordsP2);
    Calculator::calcMidPoint(coordsP1, coordsP2, midpoint);
    collapsingGrid->GetPoints()->SetPoint(p1, midpoint);
    // collapsingGrid->GetPoints()->SetPoint(p2, midpoint);
    for (auto ncell : collapse->ncells) {
        vtkIdList *ncellIds = collapsingGrid->GetCell(ncell)->GetPointIds();
        if (ncellIds->IsId(p2) != -1) {
            ncellIds->SetId(ncellIds->IsId(p2), p1);
            collapsingGrid->ReplaceCell(ncell, (int)ncellIds->GetNumberOfIds(), ncellIds->GetPointer(0));
        }
    }
    removedPointIds.insert(p2);
    vtkSmartPointer<vtkIdTypeArray> removePoints = vtkSmartPointer<vtkIdTypeArray>::New();
    removePoints->SetNumberOfComponents(1);
    removePoints->InsertNextValue(p2);
    vtkSmartPointer<vtkSelectionNode> pointNode = vtkSmartPointer<vtkSelectionNode>::New();
    pointNode->SetFieldType(vtkSelectionNode::POINT);
    pointNode->SetContentType(vtkSelectionNode::INDICES);
    pointNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);
    pointNode->SetSelectionList(removePoints);
    vtkSmartPointer<vtkSelection> pointSelection = vtkSmartPointer<vtkSelection>::New();
    pointSelection->AddNode(pointNode);
    vtkSmartPointer<vtkExtractSelection> pointExtraction = vtkSmartPointer<vtkExtractSelection>::New();
    //  ---- END NCELLS ----
    
    
    //  ---- DYNAMIC TEST ----
    //  ...
    //  Intersection Test
    //  ...
    //  ---- END DYNAMIC TEST ----
    
    
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
    }
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);
    selectionNode->SetSelectionList(icellIds);
    
    selection->AddNode(selectionNode);
    extract->SetInputData(0, collapsingGrid);
    extract->SetInputData(1, selection);
    extract->Update();
    newGrid->ShallowCopy(extract->GetOutput());
    /* ----- END ICELLS ----- */
    pointExtraction->SetInputData(0, newGrid);
    pointExtraction->SetInputData(1, pointSelection);
    pointExtraction->Update();
    collapsingGrid->ShallowCopy(pointExtraction->GetOutput());
    newGrid->SetPoints(collapsingGrid->GetPoints());
    return newGrid;
}

int main(int argc, const char * argv[]) {
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/dambreak_merged4x.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/dambreak4x.vtk";
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/wehr_clipped.vtk";
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    
    if (reader->IsFileUnstructuredGrid()) {
        
        std::cout << "input file is unstructured grid" << std::endl;
        
        // initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> simpleGrid = generateSimpleGrid();
        // std::cout << "_original_ number of cells: " << gridReader->GetOutput()->GetNumberOfCells() << std::endl;
        
        // -----------------
        // quad to tetra
        // -----------------
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        // triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->SetInputData(simpleGrid);
        triangleFilter->Update();
        std::cout << "_after triangulation_ number of cells: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();
        
        
        // --------------------
        // Initiate algorithm
        // --------------------
        std::vector<EdgeCollapse> prio_q = recalculateCosts(tetraGrid);
        // std::set<vtkIdType> removedPointIds;
        // start iteration over edge collapses
        for (int i = 0; i < 300; i++) {
            
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
                
                prio_q = recalculateCosts(tetraGrid);
                
            }
        }
        
        writeUgrid(tetraGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20151217_with-ncell.vtu");
        // startRenderingGrid(collapsingGrid);
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


