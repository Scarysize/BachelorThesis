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

std::list<Vertex> getSurfaceVertices(vtkUnstructuredGrid *grid) {
    // get the surface of the grid
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(grid);
    surfaceFilter->Update();
    
    vtkSmartPointer<vtkAppendFilter> appender = vtkSmartPointer<vtkAppendFilter>::New();
    appender->SetInputConnection(surfaceFilter->GetOutputPort());
    appender->Update();
    
    writeUgrid(appender->GetOutput(),"/Volumes/EXTERN/Bachelor Arbeit/20151207_surface-grid.vtu");
    
    std::list<Vertex> surfaceVertices;
    vtkPoints *points = surfaceFilter->GetOutput()->GetPoints();
    for (vtkIdType point = 0; point < points->GetNumberOfPoints(); point++) {
        double coords[3];
        points->GetPoint(point, coords);
        surfaceVertices.push_back(*new Vertex(coords[0], coords[1], coords[2]));
    }
    std::cout << surfaceVertices.size() << " surface vertices" << std::endl;
    return surfaceVertices;
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
    std::list<Vertex> surfaceVertices = getSurfaceVertices(tetraGrid);
    vtkSmartPointer<vtkEdgeTable> edges = vtkSmartPointer<vtkEdgeTable>::New();

    double text[3];
    tetraGrid->GetPoints()->GetPoint(66049, text);
    Vertex *test = new Vertex(text[0], text[1], text[2]);
    
    for (auto surfP : surfaceVertices) {
        if (test->hasSameCoords(surfP)) {
            std::cout << "haeoigfjweio" << std::endl;
        }
    }
    
    // Iterate over tetrahedras
    for(vtkIdType c = 0; c < 250000; c++){
        
        if ((int)c%10000 == 0) {
            std::cout << "c = " << c << std::endl;
        }
        // Iterate over cell edges
        for (vtkIdType e = 0; e < tetraGrid->GetCell(c)->GetNumberOfEdges(); e++) {
            
            vtkIdType pointA = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(0);
            vtkIdType pointB = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(1);
            
            //check if edge was already traversed
            if (edges->IsEdge(pointA, pointB) == -1) {
                
                
                /* START check for sharp vertices/edges */
                double coordsA[3];
                double coordsB[3];
                tetraGrid->GetCell(c)->GetPoints()->GetPoint(pointA, coordsA);
                tetraGrid->GetCell(c)->GetPoints()->GetPoint(pointB, coordsB);
                Vertex *vertexA = new Vertex(coordsA[0], coordsA[1], coordsA[2]);
                Vertex *vertexB = new Vertex(coordsB[0], coordsB[1], coordsB[2]);
                bool aIsSharp = false;
                bool bIsSharp = false;
                for (auto surfaceVertex : surfaceVertices) {
                    if (surfaceVertex.hasSameCoords(*vertexA)) {
                        aIsSharp = true;
                    }
                    if (surfaceVertex.hasSameCoords(*vertexB)) {
                        bIsSharp = true;
                    }
                    if (aIsSharp && bIsSharp) {
                        break;
                    }
                }
                if (aIsSharp && bIsSharp) {
                    std::cout << "surface edge" << std::endl;
                } else if (aIsSharp && !bIsSharp) {
                    std::cout << "point a on surface: " << pointA << std::endl;
                } else if (bIsSharp && !aIsSharp) {
                    std::cout << "point b on surface: "<< pointB << std::endl;
                } else {
                    std::cout << "is interior edge" << std::endl;
                }
                /* END check for sharp vertices */
                if ((aIsSharp && bIsSharp) || (!aIsSharp && !bIsSharp)) {
                    std::set<vtkIdType> icells = CostCalculator::getIntroducedTetras((int) e, c, tetraGrid);
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
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
    
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
        
        for (int i = 0; i < 50; i++) {
            if (!prio_q.empty()) {
                EdgeCollapse top = prio_q.top();
                doCollapse(&top, collapsingGrid);
                //prio_q = recalculateCosts(collapsingGrid);
                prio_q.pop();
                std::cout << "_after first simplification_: " << collapsingGrid->GetNumberOfCells() << std::endl;
            }
        }
        
        
        
        writeUgrid(collapsingGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20151207_with-ncell.vtu");
        // startRenderingGrid(collapsingGrid);
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


