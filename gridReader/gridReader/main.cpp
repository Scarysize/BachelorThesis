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

bool isInnerVertex(double solidAngle) {
    if (solidAngle >= 3.1) {
        return true;
    }
    return false;
}

bool isCorner(double solidAngle) {
    if (solidAngle <= M_PI) {
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
            vertexIds.insert(seed);
        } else if (!isInnerVertex(angle)) {
            vertexIds.insert(seed);
        }
    }
    
    return vertexIds;
}


/*!
 Does the cost calucation step and builds up the edge collapse heap in the process
 \param tetraGrid The unstructured grid to calculate the cost and edged collapses from
 \returns The priority queue containing all possible edge collapses
 */
std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, CompareCost> recalculateCosts(vtkUnstructuredGrid *tetraGrid) {
    std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, CompareCost> prio_q;
    // ----------------------
    // EXTRACT FEATURE EDGES
    // ----------------------
    std::set<vtkIdType> surfaceIds = getSurfaceVertexIds(tetraGrid, getSurfaceVertices(tetraGrid));
    vtkSmartPointer<vtkEdgeTable> edges = vtkSmartPointer<vtkEdgeTable>::New();
    
    // Iterate over tetrahedras
    for(vtkIdType c = 0; c < 5000; c++){
        
        if ((int)c%1000 == 0) {
            std::cout << "c = " << c << std::endl;
        }
        // Iterate over cell edges
        for (vtkIdType e = 0; e < tetraGrid->GetCell(c)->GetNumberOfEdges(); e++) {
            
            vtkIdType pointA = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(0);
            vtkIdType pointB = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(1);
            
            //check if edge was already traversed
            if (edges->IsEdge(pointA, pointB) == -1) {
                
                std::set<vtkIdType> icells = CostCalculator::getIntroducedTetras((int) e, c, tetraGrid);
                
                /* START check for sharp vertices/edges */
                bool aIsSharp = false;
                bool bIsSharp = false;
                bool hasIcellProblem = false;
                
                if (surfaceIds.find(pointA) != surfaceIds.end()) {
                    //std::cout << "point is on surface(a): " << pointA << std::endl;
                    aIsSharp = true;
                }
                if (surfaceIds.find(pointB) != surfaceIds.end()) {
                    //std::cout << "point is on surface(b): " << pointB << std::endl;
                    bIsSharp = true;
                }
                
                for (auto icell : icells) {
                    if (surfaceIds.find(icell) != surfaceIds.end()) {
                        hasIcellProblem = true;
                    }
                }

                if (((aIsSharp && bIsSharp) || (!aIsSharp && !bIsSharp)) && !hasIcellProblem) {
                    std::set<vtkIdType> ncells = CostCalculator::getNonVanishingTetras((int) e, c, tetraGrid);
                    double volumeCost = CostCalculator::calcVolumeCost((int) e, c, 500, &ncells, &icells , tetraGrid);
                    double scalarCost = CostCalculator::calcScalarCost(tetraGrid->GetCell(c)->GetEdge((int)e), 1, tetraGrid);
                    double collapseCost = volumeCost + scalarCost;
                    EdgeCollapse collapse = *new EdgeCollapse((int)e, c, collapseCost, icells, ncells);
                    edges->InsertEdge(pointA, pointB);
                    prio_q.push(collapse);
                }
            }
            
        }
        
    }
    
    return prio_q;
}

/*!
 Does one edge collapse. Changes the coords of ncells accordingly and removes icells from the grid.
 \param collapse The EdgeCollapse to be executed.
 \param collapsingGrid The unstructured grid the collapse should be executed on.
 \returns The unstructured grid after the collapse.
 */
void doCollapse(EdgeCollapse *collapse, vtkUnstructuredGrid *collapsingGrid) {
    
    /* manipulate ncells */
    vtkIdType p1 = collapsingGrid->GetCell(collapse->getTetraId())->GetEdge(collapse->getEdgeId())->GetPointId(0);
    vtkIdType p2 = collapsingGrid->GetCell(collapse->getTetraId())->GetEdge(collapse->getEdgeId())->GetPointId(1);
    
    double midpoint[3];
    collapse->calcCollapsePoint(collapsingGrid, midpoint);
    collapsingGrid->GetPoints()->SetPoint(p1, midpoint);
    collapsingGrid->GetPoints()->SetPoint(p2, midpoint);
    
    vtkSmartPointer<vtkPoints> oldPoints = collapsingGrid->GetPoints();
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    
    for (vtkIdType pointId = 0; pointId < oldPoints->GetNumberOfPoints(); pointId++) {
        double point[3];
        oldPoints->GetPoint(pointId, point);
        newPoints->InsertNextPoint(point);
    }
    
    /* remove icells */
    vtkSmartPointer<vtkIdTypeArray> removeCells = vtkSmartPointer<vtkIdTypeArray>::New();
    removeCells->SetNumberOfComponents(1);
    std::set<vtkIdType> removeCellIds = collapse->icells;
    for (auto id : removeCellIds) {
        removeCells->InsertNextValue(id);
    }
    
    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1); //invert selection
    selectionNode->SetSelectionList(removeCells);
    
    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);
    
    vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, collapsingGrid);
    extractSelection->SetInputData(1, selection);
    
    extractSelection->Update();
    collapsingGrid->ShallowCopy(extractSelection->GetOutput());
}

int main(int argc, const char * argv[]) {
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    if (reader->IsFileUnstructuredGrid()) {
        
        reader = NULL;
        
        std::cout << "input file is unstructured grid" << std::endl;
        
        //Initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        std::cout << "_original_ number of cells: " << gridReader->GetOutput()->GetNumberOfCells() << std::endl;
        
        
        // -----------------
        // quad to tetra
        // -----------------
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->Update();
        std::cout << "_after triangulation_ number of cells: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();
        gridReader = NULL; //don´t need this anymore
        
        std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, CompareCost> prio_q = recalculateCosts(tetraGrid);
        
        vtkSmartPointer<vtkUnstructuredGrid> collapsingGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        collapsingGrid->ShallowCopy(triangleFilter->GetOutput());
        
        std::cout << "prio queue: " << prio_q.size() << std::endl;
        
        for (int i = 0; i < 100; i++) {
            if (!prio_q.empty()) {
                EdgeCollapse top = prio_q.top();
                doCollapse(&top, collapsingGrid);
                //prio_q = recalculateCosts(collapsingGrid);
                prio_q.pop();
            }
        }
        std::cout << "_after first simplification_: " << collapsingGrid->GetNumberOfCells() << std::endl;
        
        
        
        writeUgrid(collapsingGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20151207_with-ncell.vtu");
        // startRenderingGrid(collapsingGrid);
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


