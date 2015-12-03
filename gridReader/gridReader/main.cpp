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


#include <vtkIdList.h>
#include <vtkIdTypeArray.h>

#include <vtkAlgorithmOutput.h>

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

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
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "Calculator.h"
#include "CostCalculator.hpp"
#include "EdgeCollapse.hpp"
#include "Helper.hpp"


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

void startRendering(vtkAlgorithmOutput* grid) {
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

int main(int argc, const char * argv[]) {
    vtkIndent *indent = vtkIndent::New(); //for PrintSelf calls
    
    
    
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    // std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Wehr/01_inter_wehr/VTK/01_inter_wehr_LES_SpalartAllmarasDDES_12891.vtk";
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    if (reader->IsFileUnstructuredGrid()) {
        
        struct CompareCost {
            bool operator()(EdgeCollapse &col1, EdgeCollapse &col2) {
                return col1.getCost() > col2.getCost();
            }
        };
        reader = NULL;
        
        std::cout << "input file is unstructured grid" << std::endl;
        
        std::priority_queue<EdgeCollapse, std::vector<EdgeCollapse>, CompareCost> prio_q;
        
        //Initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        // --------------
        // ORIGINAL GRID
        // --------------
        // vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // ugrid = gridReader->GetOutput();
        std::cout << "_original_ number of cells: " << gridReader->GetOutput()->GetNumberOfCells() << std::endl;
        // free(ugrid);
        
        // -----------------
        // TETRAHEDRAL GRID
        // -----------------
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->Update();
        std::cout << "_after triangulation_ number of cells: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        
        vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
        geoFilter->SetInputConnection(triangleFilter->GetOutputPort());
        geoFilter->Update();
        vtkSmartPointer<vtkFeatureEdges> features = vtkSmartPointer<vtkFeatureEdges>::New();
        features->SetInputConnection(geoFilter->GetOutputPort());
        features->BoundaryEdgesOn();
        features->FeatureEdgesOn();
        features->ManifoldEdgesOn();
        features->NonManifoldEdgesOn();
        features->Update();
        
        vtkSmartPointer<vtkEdgeTable> featureTable = vtkSmartPointer<vtkEdgeTable>::New();
        vtkPolyData *featureEdges = features->GetOutput();
        for (vtkIdType fe = 0; fe < featureEdges->GetNumberOfCells(); fe++) {
            vtkIdType pointA = featureEdges->GetCell(fe)->GetPointId(0);
            vtkIdType pointB = featureEdges->GetCell(fe)->GetPointId(1);
            if (featureTable->IsEdge(pointA, pointB) == -1) {
                featureTable->InsertEdge(pointA, pointB);
            }
        }
        
        std::cout << featureTable->GetNumberOfEdges() << std::endl;

        
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();
        gridReader = NULL;
        
        vtkSmartPointer<vtkEdgeTable> edges = vtkSmartPointer<vtkEdgeTable>::New();
        
        // Iterate over tetrahedras
        for(vtkIdType c = 0; c < tetraGrid->GetNumberOfCells(); c++){
            
            for (vtkIdType e = 0; e < tetraGrid->GetCell(c)->GetNumberOfEdges(); e++) {

                vtkIdType pointA = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(0);
                vtkIdType pointB = tetraGrid->GetCell(c)->GetEdge((int) e)->GetPointId(1);
                // surfacePoints->GetPoin;
                
                //check if edge was already traversed
                if (edges->IsEdge(pointA, pointB) == -1 && featureTable->IsEdge(pointA, pointB) != -1) {
                    edges->InsertEdge(tetraGrid->GetCell(c)->GetEdge((int)e)->GetPointId(0), tetraGrid->GetCell(c)->GetEdge((int)e)->GetPointId(1));
                    std::set<vtkIdType> icells = CostCalculator::getIntroducedTetras((int) e, c, tetraGrid);
                    std::set<vtkIdType> ncells = CostCalculator::getNonVanishingTetras((int) e, c, tetraGrid);
                    double volumeCost = CostCalculator::calcVolumeCost((int) e, c, 500, &ncells, &icells , tetraGrid);
                    double scalarCost = CostCalculator::calcScalarCost(tetraGrid->GetCell(c)->GetEdge((int)e), 1, tetraGrid);
                    double collapseCost = volumeCost + scalarCost;
                    EdgeCollapse *collapse = new EdgeCollapse((int)e, c, collapseCost, &icells, &ncells);
                    prio_q.push(*collapse);
                }
                
            }
            
        }
        std::cout << "number of edges: " << edges->GetNumberOfEdges() << std::endl;
        
        std::cout << "size of prio_q: " << prio_q.size() << std::endl;
        
        vtkSmartPointer<vtkPoints> points = tetraGrid->GetPoints();
       for (int i = 0; i < prio_q.size(); i++) {
            EdgeCollapse collapse = prio_q.top();
            // std::cout << "cost: " << collapse.getCost() << std::endl;
            prio_q.pop();
            
            vtkIdType p1 = tetraGrid->GetCell(collapse.getTetraId())->GetEdge(collapse.getEdgeId())->GetPointId(0);
            vtkIdType p2 = tetraGrid->GetCell(collapse.getTetraId())->GetEdge(collapse.getEdgeId())->GetPointId(1);
            
            double midpoint[3];
            collapse.calcCollapsePoint(tetraGrid, midpoint);
            points->SetPoint(p1, midpoint);
            points->SetPoint(p2, midpoint);
            
        }
        
        std::cout << "_after first simplification_: " << tetraGrid->GetNumberOfCells() << std::endl;
        
        // writeUgrid(tetraGrid, "/Volumes/EXTERN/Bachelor Arbeit/test_20151203.vtu");
        // startRendering(triangleFilter->GetOutputPort());
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


