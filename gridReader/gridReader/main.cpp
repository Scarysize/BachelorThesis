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
#include <iostream>
#include <vtkSmartPointer.h>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridGeometryFilter.h>

#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>

#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkDelaunay3D.h>

#include <vtkGeometryFilter.h>

#include <vtkQuadricDecimation.h>
#include <vtkDecimatePro.h>

#include <vtkFieldDataToAttributeDataFilter.h>

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>


vtkSmartPointer<vtkActor> delaunay3d(std::string filename) {
    //Read .vtk file
    vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    gridReader->SetFileName(filename.c_str());
    gridReader->Update();
    
    vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geoFilter->SetInputData(gridReader->GetOutput());
    geoFilter->Update();
    vtkPolyData *gridToPoly = geoFilter->GetOutput();
    
    //Get clean poly data --> removes duplicates
    vtkSmartPointer<vtkCleanPolyData> cleanpoly = vtkSmartPointer<vtkCleanPolyData>::New();
    //cleanpoly->SetInputData(gridToPoly);
    
    //Generate 3d triangle mesh
    vtkSmartPointer<vtkDelaunay3D> delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay3D->SetInputData(gridToPoly);
    
    //WriteDataSet(delaunay3D->GetOutput(), "/Volumes/EXTERN/Bachelor Arbeit/delaunay3d.vtp");
    
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetInputConnection(delaunay3D->GetOutputPort());
    writer->SetFileName("/Volumes/EXTERN/delaunay3d.vtu");
    writer->Write();
    
    vtkSmartPointer<vtkDataSetMapper> delaunayMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    delaunayMapper->SetInputConnection(delaunay3D->GetOutputPort());
    
    vtkSmartPointer<vtkActor> delaunayActor = vtkSmartPointer<vtkActor>::New();
    delaunayActor->SetMapper(delaunayMapper);
    delaunayActor->GetProperty()->SetColor(1,0,0);
    
    return delaunayActor;
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
        
        vtkSmartPointer<vtkActor> delaunay = delaunay3d(inputFilename);
        
        //Initialize UnstructuredGridReader
        vtkSmartPointer<vtkUnstructuredGridReader> gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        gridReader->SetFileName(inputFilename.c_str());
        gridReader->Update();
        
        vtkSmartPointer<vtkUnstructuredGrid> ugrid = gridReader->GetOutput();
        std::cout << "cell type 0; " << ugrid->GetCellType(0) << std::endl;
        std::cout << "cell type 100; " << ugrid->GetCellType(100) << std::endl;
        std::cout << "cell type 1000; " << ugrid->GetCellType(1000) << std::endl;
        std::cout << "cell type 6000; " << ugrid->GetCellType(6000) << std::endl;
        
        
        //Convert grid to PolyData for decimation
        /*vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
        geoFilter->SetInputData(ugrid);
        geoFilter->Update();
        vtkPolyData *poly = geoFilter->GetOutput();
        
        std::cout << "polydata cells:  " << poly->GetNumberOfCells() << std::endl;
        std::cout << "polydata points: " << poly->GetNumberOfPoints() << std::endl;
        
        //Decimation
        vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
        decimate->SetInputData(poly);
        decimate->SetTargetReduction(.1);
        decimate->Update();
        
        vtkSmartPointer<vtkPolyData> decimated = vtkSmartPointer<vtkPolyData>::New();
        decimated->ShallowCopy(decimate->GetOutput());
        
        std::cout << "decimated number of cells: " << decimated->GetNumberOfCells() << std::endl;*/
        
        
        //Field data --> to Cells
        vtkSmartPointer<vtkFieldDataToAttributeDataFilter> toCellFilter = vtkSmartPointer<vtkFieldDataToAttributeDataFilter>::New();
        toCellFilter->SetInputConnection(gridReader->GetOutputPort());
        toCellFilter->SetInputFieldToCellDataField();
        toCellFilter->SetOutputAttributeDataToCellData();
        
        
        toCellFilter->SetScalarComponent(0, "alpha.water", 0);

        
        //Mapper
        vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputConnection(toCellFilter->GetOutputPort());
        
        //Actor
        vtkSmartPointer<vtkActor> ugridActor = vtkSmartPointer<vtkActor>::New();
        ugridActor->SetMapper(mapper);
        
        
        //Renderer
        vtkSmartPointer<vtkRenderer> ugridRenderer = vtkSmartPointer<vtkRenderer>::New();
        ugridRenderer->AddActor(delaunay);
        ugridRenderer->SetBackground(0.1, 0.2, 0.4);
        
        //Render Window & Interactor
        vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
        renWin->AddRenderer(ugridRenderer);
        renWin->SetSize(1024, 1024);
        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(renWin);
        iren->Initialize();
        iren->Start();
        
        //std::cout << "number of cells: " << ugrid->GetNumberOfCells() << std::endl;
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "output is structured grid" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


