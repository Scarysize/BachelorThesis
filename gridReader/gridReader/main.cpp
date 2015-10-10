//
//  main.cpp
//  gridReader
//
//  Created by Franz Neubert on 06/10/15.
//  Copyright © 2015 Franz Neubert. All rights reserved.
//

#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL);
#include <vtkGenericDataObjectReader.h>
#include <iostream>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

void visualize(vtkSmartPointer<vtkPolyDataMapper> mapper) {
    
 
}

int main(int argc, const char * argv[]) {
    
    
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/XCode/gridReader/gridReader/damnbreak_708.vtk";
    //std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/01_inter_RHG_BHQ1_SA_mesh01_0.vtk";
    std::string inputFilename = "/Volumes/EXTERN/Bachelor Arbeit/OpenFOAM_Daten/Rinne/inter_RHG/VTK/wall_Rinne/wall_Rinne_0.vtk";
    
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    
    if(reader->IsFilePolyData())
    {
        std::cout << "output is polydata" << std::endl;
        vtkPolyData *ouput = reader->GetPolyDataOutput();
        
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
        
        std::cout << "output has " << ouput->GetNumberOfPoints() << " points." << std::endl;
    } else if (reader->IsFileUnstructuredGrid()) {
        std::cout << "ouput is unstructured grid" << std::endl;
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->Update();
        
        vtkSmartPointer<vtkUnstructuredGridGeometryFilter> geometryFilter = vtkSmartPointer<vtkUnstructuredGridGeometryFilter>::New();
        geometryFilter->SetInputConnection(reader->GetOutputPort());
        geometryFilter->Update();
        
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(geometryFilter->GetOutputPort());
        
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
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "output is structured grid" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


