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

#include <iostream>
#include <vtkSmartPointer.h>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridGeometryFilter.h>

#include <vtkHexahedron.h>

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

#include <vtkQuadricDecimation.h>
#include <vtkDecimatePro.h>

#include <vtkFieldDataToAttributeDataFilter.h>

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

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

void startRendering(vtkAlgorithmOutput *out) {
    //Mapper
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(out);
    
    //Actor
    vtkSmartPointer<vtkActor> ugridActor = vtkSmartPointer<vtkActor>::New();
    ugridActor->SetMapper(mapper);
    
    
    //Renderer
    vtkSmartPointer<vtkRenderer> ugridRenderer = vtkSmartPointer<vtkRenderer>::New();
    ugridRenderer->AddActor(ugridActor);
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
    
    
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Dambreak/00_damBreak_2d/01_inter/VTK/01_inter_50.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/wall_Rinne/wall_Rinne_0.vtk";
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    if(reader->IsFilePolyData())
    {
        std::cout << "output is polydata" << std::endl;
        //vtkPolyData *ouput = reader->GetPolyDataOutput();
        
        //MAPPER
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());
        
        //ACTOR
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        
        //RENDERER
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(renderWindow);
        
        renderer->AddActor(actor);
        renderer->SetBackground(0.1, 0.2, 0.4);
        renderWindow->SetSize(1024, 1024);
        renderWindow->Render();
        
        iren->Start();
        
        //std::cout << "output has " << ouput->GetNumberOfPoints() << " points." << std::endl;
    } else if (reader->IsFileUnstructuredGrid()) {
        std::cout << "ouput is unstructured grid" << std::endl;
        
        vtkIndent *indent = vtkIndent::New(); //for PrintSelf calls
        
        //Initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> ugrid = gridReader->GetOutput();
        std::cout << ugrid->GetNumberOfCells() << std::endl;
        
        vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        triangleFilter->SetInputConnection(gridReader->GetOutputPort());
        triangleFilter->Update();
        
        std::cout << triangleFilter->GetOutput()->GetNumberOfCells() << std::endl;
        writeUgrid(triangleFilter->GetOutput(), "/Volumes/EXTERN/trifilter_dambreak.vtu");
        
        startRendering(triangleFilter->GetOutputPort());
        
        /*vtkSmartPointer<vtkPolyDataMapper> pMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        pMapper->SetInputData(gridToPoly);
        
        vtkSmartPointer<vtkActor> pActor = vtkSmartPointer<vtkActor>::New();
        pActor->SetMapper(pMapper);
        
        vtkSmartPointer<vtkRenderer> pRen = vtkSmartPointer<vtkRenderer>::New();
        pRen->AddActor(pActor);
        
        vtkSmartPointer<vtkRenderWindow> pWin = vtkSmartPointer<vtkRenderWindow>::New();
        pWin->AddRenderer(pRen);
        pWin->SetSize(1024, 1024);
        vtkSmartPointer<vtkRenderWindowInteractor> pIn = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        pIn->SetRenderWindow(pWin);
        pIn->Initialize();
        pIn->Start();*/
        
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "output is structured grid" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


