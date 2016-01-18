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

#include <vtkTetra.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkCell.h>
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

#include <vtkGeometryFilter.h>
#include <vtkAppendFilter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkExtractEdges.h>
#include <vtkFeatureEdges.h>
#include <vtkExtractUnstructuredGrid.h>

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
#include "DynamicTest.hpp"
#include "Cell.hpp"
#include "Connectivity.hpp"
#include "GridReducer.hpp"
#include "CostCalculations.hpp"

struct CompareCell {
    bool operator()(Cell &cell1, Cell &cell2) {
        return cell1.id == cell2.id;
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

void startRenderingGrid(vtkUnstructuredGrid *grid, vtkRenderWindow *renWin) {
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
    ugridActor->GetProperty()->SetOpacity(0.9);
    
    
    //Renderer
    vtkSmartPointer<vtkRenderer> ugridRenderer = vtkSmartPointer<vtkRenderer>::New();
    // ugridRenderer->AddActor(featureActor);
    ugridRenderer->AddActor(ugridActor);
    ugridRenderer->AddActor(edgeActor);
    ugridRenderer->SetBackground(0.1, 0.2, 0.4);
    
    //Render Window & Interactor
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

bool isInner(double solidAngle) {
    if (solidAngle >= 12.566370614359176) {
        return true;
    }
    return false;
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



void iterateTetras(vtkUnstructuredGrid *tetraGrid,
                   std::vector<EdgeCollapse> *prio_q,
                   std::set<vtkIdType> surfacePoints,
                   std::set<vtkIdType> featurePoints) {
}

/*!
 Does the cost calucation step and builds up the edge collapse heap in the process
 \param tetraGrid The unstructured grid to calculate the cost and edged collapses from
 \returns The priority queue containing all possible edge collapses
 */


/*!
 Does one edge collapse. Changes the coords of ncells accordingly and removes icells from the grid.
 \param collapse The EdgeCollapse to be executed.
 \param collapsingGrid The unstructured grid the collapse should be executed on.
 \returns The unstructured grid after the collapse.
 */


int main(int argc, const char * argv[]) {
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/dambreak_merged4x.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/dambreak4x.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/wehr_clipped.vtk";
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    
    if (reader->IsFileUnstructuredGrid()) {
        
        std::cout << "INFO: input file is unstructured grid" << std::endl;
        
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
        triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        // triangleFilter->SetInputData(simpleGrid);
        triangleFilter->Update();
        std::cout << "INFO: cells after triangulation: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();
        
        std::vector<Cell*> cells = Cell::cellsFromGrid(tetraGrid);
        std::vector<Vertex*> vertices = Vertex::verticesFromGrid(tetraGrid);
        std::cout << "DONE: data retrieval" << std::endl;
        
        GridReducer *reducer = new GridReducer(cells, vertices);
        reducer->run(&CostCalculations::calcEdgeLengthCost);
        
        std::vector<Cell*> postColCells = reducer->getCells();
        std::vector<Vertex*> postColVertices = reducer->getVertices();
        vtkSmartPointer<vtkUnstructuredGrid> postColGrid = Helper::makeGrid(postColCells, postColVertices);
        writeUgrid(postColGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20160118.vtu");
        
        std::cout << "INFO: done -----------------" << std::endl;
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


