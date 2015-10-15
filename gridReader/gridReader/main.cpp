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
#include <vtkCell.h>
#include <vtkPointData.h>

#include <vtkDoubleArray.h>

#include <vtkAlgorithmOutput.h>

#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkDelaunay3D.h>

#include <vtkGeometryFilter.h>

#include <vtkAppendFilter.h>

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

vtkPolyData* ugridToPolyData (vtkUnstructuredGrid *grid) {
    vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geoFilter->SetInputData(grid);
    geoFilter->Update();
    return geoFilter->GetOutput();
}

vtkUnstructuredGrid* datasetToUgrid (vtkDataSet *dataset) {
    //Append filter: Dataset --> to UnstructuredGrid
    vtkSmartPointer<vtkAppendFilter> appender = vtkSmartPointer<vtkAppendFilter>::New();
    appender->AddInputData(dataset);
    appender->Update();
    
    return appender->GetOutput();
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
        
        //Field data --> to Cells
        vtkSmartPointer<vtkFieldDataToAttributeDataFilter> toCellFilter = vtkSmartPointer<vtkFieldDataToAttributeDataFilter>::New();
        toCellFilter->SetInputConnection(gridReader->GetOutputPort());
        toCellFilter->SetInputFieldToCellDataField();
        toCellFilter->SetOutputAttributeDataToCellData();
        
        /*vtkDataArray *arr = ugrid->GetPointData()->GetArray("alpha.water");
        std::cout << arr->GetSize() << std::endl;*/
        
        toCellFilter->SetScalarComponent(0, "alpha.water", 0);
        //toCellFilter->SetScalarComponent(0, "U", 0);
        
        startRendering(toCellFilter->GetOutputPort());
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "output is structured grid" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


