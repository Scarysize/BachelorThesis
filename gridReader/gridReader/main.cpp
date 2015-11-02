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

#include <vtkGenericDataObjectReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include <vtkIdList.h>
#include <vtkIdTypeArray.h>

#include <vtkAlgorithmOutput.h>

#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkDelaunay3D.h>

#include <vtkGeometryFilter.h>

#include <vtkAppendFilter.h>

#include <vtkDataSetTriangleFilter.h>

#include <vtkExtractEdges.h>

#include <vtkQuadricDecimation.h>
#include <vtkDecimatePro.h>

#include <vtkFieldDataToAttributeDataFilter.h>

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "Calculator.h"

void printVector(double x[3]) {
    std::vector<double> vector(x, x + sizeof(double[3]) / sizeof(double));
    for(std::vector<double>::const_iterator i = vector.begin(); i != vector.end(); ++i) {
        std::cout << *i << ' ';
    }
    std::cout << std::endl;
}

double calcCellGradient(vtkTetra *cell, vtkUnstructuredGrid *tetraGrid) {
    double sum;
    double gradientSum = 0;
    
    // Get Point Coordinates
    vtkPoints* points = cell->GetPoints();
    double A[3];
    double B[3];
    double C[3];
    double D[3];
    points->GetPoint(0, A);
    points->GetPoint(1, B);
    points->GetPoint(2, C);
    points->GetPoint(4, D);
    
    // Calculate necessary vectors for gradients on each vertex
    double BD[3];
    double BC[3];
    double AC[3];
    double AD[3];
    double AB[3];
    Calculator::subtractVectors(B, D, BD);
    Calculator::subtractVectors(B, C, BC);
    Calculator::subtractVectors(A, C, AC);
    Calculator::subtractVectors(A, D, AD);
    Calculator::subtractVectors(A, B, AB);
    
    // Caluclate the gradients on each vertex (via cross product)
    double gradA[3];
    double gradB[3];
    double gradC[3];
    double gradD[3];
    Calculator::calcCrossProduct(BD, BC, gradA);
    Calculator::calcCrossProduct(AC, AD, gradB);
    Calculator::calcCrossProduct(AD, AB, gradC);
    Calculator::calcCrossProduct(AC, AB, gradD);
    
    
    gradientSum += Calculator::calcVectorLength(gradA);
    gradientSum += Calculator::calcVectorLength(gradB);
    gradientSum += Calculator::calcVectorLength(gradC);
    gradientSum += Calculator::calcVectorLength(gradD);
    
    std::cout << "gradient sum: " << gradientSum << std::endl;
    
    sum = (gradientSum * fabs(cell->ComputeVolume(A, B, C, D)))/4;
    
    std::cout << "unweighted gradient: " << sum << std::endl;
    
    return sum;
}



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
 Checks if two sets of cell neighbours have a common neighbour.
 */
bool shareNeighbours(std::list<vtkIdType> currentCellNeighbourIds, std::list<vtkIdType> lastMergedCellIds) {
    for(std::list<vtkIdType>::iterator currentIt = currentCellNeighbourIds.begin(); currentIt != currentCellNeighbourIds.end(); currentIt++) {
        for (std::list<vtkIdType>::iterator lastMergedIt = lastMergedCellIds.begin(); lastMergedIt != lastMergedCellIds.end(); lastMergedIt++) {
            if (*currentIt == *lastMergedIt) {
                return true;
            }
        }
    }
    return false;
}

/*
 Returns a id list of all neighbour cells of a given cell in a given grid.
 */
void getNeighbourCellIds (vtkUnstructuredGrid *ugrid, vtkIdType cellId, std::list<vtkIdType> *neighbours){
    
    vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
    ugrid->GetCellPoints(cellId, cellPointIds);
    
    for(vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
        //Check for every cell point: which cells (in the grid) share this point --> neighbouring cells
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        idList->InsertNextId(cellPointIds->GetId(i));
        vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
        ugrid->GetCellNeighbors(cellId, idList, neighbourCellIds);
        //save the ids of the neighbour cells, ommit duplicates
        for(vtkIdType j = 0; j < neighbourCellIds->GetNumberOfIds(); j++) {
            bool inArray = (std::find(neighbours->begin(), neighbours->end(), neighbourCellIds->GetId(j)) != neighbours->end());
            if(!inArray) {
                neighbours->push_back(neighbourCellIds->GetId(j));
            }
        }
    }
    
    
}

void printCellNeighbours (vtkIdType cellId, std::list<vtkIdType> *neighbours) {
    std::cout << "Neighbours current cell " << cellId << std::endl;
    for(std::list<vtkIdType>::iterator it1 = neighbours->begin(); it1 != neighbours->end(); it1++) {
        std::cout << " " << *it1;
    }
    std::cout << std::endl << "- - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
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
    ugridRenderer->AddActor(ugridActor);
    ugridRenderer->AddActor(edgeActor);
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
    
    
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/wall_Rinne/wall_Rinne_0.vtk";
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    if (reader->IsFileUnstructuredGrid()) {
        std::cout << "input file is unstructured grid" << std::endl;
        
        
        //Initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        // --------------
        // ORIGINAL GRID
        // --------------
        vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        ugrid = gridReader->GetOutput();
        std::cout << "_original_ number of cells: " << ugrid->GetNumberOfCells() << std::endl;
        
        // -----------------
        // TETRAHEDRAL GRID
        // -----------------
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->Update();
        std::cout << "_after triangulation_ number of cells: " << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        
        vtkSmartPointer<vtkUnstructuredGrid>  tetraGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        tetraGrid = triangleFilter->GetOutput();
        
        // Iterate over tetrahedras
        vtkIdType cellId;
        for(vtkIdType c = 0; c < 50 /* tetraGrid->GetNumberOfCells() */; c++){
            cellId = c;
            vtkTetra *tetra = vtkTetra::SafeDownCast(tetraGrid->GetCell(cellId));
            calcCellGradient(tetra, tetraGrid);
            //std::cout << "cell type: " << tetraGrid->GetCell(cellId)->GetCellType() << std::endl;
        }
        //startRendering(triangleFilter->GetOutputPort());
        
        
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "input file is structured grid" << std::endl;
    } else {
        std::cerr << "error reading file" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


