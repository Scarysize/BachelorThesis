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
        
        //Field data --> to Cells
        vtkSmartPointer<vtkFieldDataToAttributeDataFilter> toCellFilter = vtkSmartPointer<vtkFieldDataToAttributeDataFilter>::New();
        toCellFilter->SetInputConnection(gridReader->GetOutputPort());
        toCellFilter->SetInputFieldToCellDataField();
        toCellFilter->SetOutputAttributeDataToCellData();
        
        /*vtkDataArray *arr = ugrid->GetPointData()->GetArray("alpha.water");
         vtkFloatArray *floatArr = vtkFloatArray::SafeDownCast(arr);
         if(floatArr){
         std::cout << "is array" << std::endl;
         for(int i = 0; i < floatArr->GetSize(); i++) {
         float alphaWater;
         alphaWater = floatArr->GetValue(i);
         if(i%100 == 0){
         std::cout << alphaWater << std::endl;
         }
         }
         }*/
        
        //Iterator over all grid cells
        vtkIdType cellId;
        vtkIdType lastCellId = -1;
        int mergeCounter = 0;
        std::list<vtkIdType> lastCellNeighbours;
        vtkSmartPointer<vtkPoints> mergePoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> hexs = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkUnstructuredGrid> mergeGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        
        for(vtkIdType c = 0; c < ugrid->GetNumberOfCells(); c++){
            cellId = c;

            std::list<vtkIdType> neighbours;
            getNeighbourCellIds(ugrid, cellId, &neighbours);
            if (shareNeighbours(neighbours, lastCellNeighbours)) {
                std::cout << cellId << " and " << lastCellId << " share neighbours" << std::endl;
            } else {
                if (neighbours.size() == 8) {
                    std::cout << "can merge: " << cellId << std::endl;
                }
            }
            
            /*
             Special case for dambreak:
             
             - INNER cells have 8 neighbours
             - BOUNDARY cells have 5 neighbours
             - CORNER cells have 3 neighbours
             
             this is due to the "2D" layout (only 1 row of cells on the z-axis
             
             ==> Only merge if cell is INNER cell
             
             TODO: Adapt this for "real" 3D Grid
            */
            if(neighbours.size() == 8 && !shareNeighbours(neighbours, lastCellNeighbours)) {
                mergeCounter++;
                lastCellId = cellId;
                lastCellNeighbours = neighbours;
                /*
                    Build a collection of all points of every neighbour cell
                */
                vtkSmartPointer<vtkPoints> neighbourhoodPoints = vtkSmartPointer<vtkPoints>::New();
                for(std::list<vtkIdType>::iterator it1 = neighbours.begin(); it1 != neighbours.end(); it1++) {
                    vtkPoints *cellPoints = ugrid->GetCell(*it1)->GetPoints();
                    for (vtkIdType p = 0; p < cellPoints->GetNumberOfPoints(); p++) {
                        double pointComponents[3];
                        cellPoints->GetPoint(p, pointComponents);
                        neighbourhoodPoints->InsertNextPoint(pointComponents);
                    }
                }
                
                /*
                    Get the outer bounds of the collection of points
                */
                double bounds[6];
                neighbourhoodPoints->ComputeBounds();
                neighbourhoodPoints->GetBounds(bounds);
                double xmin = bounds[0];
                double xmax = bounds[1];
                double ymin = bounds[2];
                double ymax = bounds[3];
                double zmin = bounds[4];
                double zmax = bounds[5];
                
                /*
                    Build a collection of corner points for the outer bounds
                */
                double p0[3] = {xmin, ymin, zmin};
                double p1[3] = {xmax, ymin, zmin};
                double p2[3] = {xmax, ymax, zmin};
                double p3[3] = {xmin, ymax, zmin};
                double p4[3] = {xmin, ymin, zmax};
                double p5[3] = {xmax, ymin, zmax};
                double p6[3] = {xmax, ymax, zmax};
                double p7[3] = {xmin, ymax, zmax};

                
                /*
                    Create a hexahedron from the new corner points
                */
                vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
                hex->GetPointIds()->SetId(0, mergePoints->InsertNextPoint(p0));
                hex->GetPointIds()->SetId(1, mergePoints->InsertNextPoint(p1));
                hex->GetPointIds()->SetId(2, mergePoints->InsertNextPoint(p2));
                hex->GetPointIds()->SetId(3, mergePoints->InsertNextPoint(p3));
                hex->GetPointIds()->SetId(4, mergePoints->InsertNextPoint(p4));
                hex->GetPointIds()->SetId(5, mergePoints->InsertNextPoint(p5));
                hex->GetPointIds()->SetId(6, mergePoints->InsertNextPoint(p6));
                hex->GetPointIds()->SetId(7, mergePoints->InsertNextPoint(p7));

                hexs->InsertNextCell(hex);
                
                mergeGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
            }
            
        }
        std::cout << "merge counter: " << mergeCounter << std::endl;
        mergeGrid->SetPoints(mergePoints);
        
        writeUgrid(mergeGrid, "/Volumes/EXTERN/Bachelor Arbeit/XCode/BachelorThesis/gridReader/mergeGrid.vtu");
        mergeGrid->PrintSelf(std::cout, *indent);
        
        vtkSmartPointer<vtkDataSetMapper> mp = vtkSmartPointer<vtkDataSetMapper>::New();
        mp->SetInputData(mergeGrid);
        vtkSmartPointer<vtkActor> ac = vtkSmartPointer<vtkActor>::New();
        ac->SetMapper(mp);
        vtkSmartPointer<vtkRenderer> rn = vtkSmartPointer<vtkRenderer>::New();
        rn->AddActor(ac);
        rn->SetBackground(.2, .3, .4);
        vtkSmartPointer<vtkRenderWindow> rw = vtkSmartPointer<vtkRenderWindow>::New();
        rw->AddRenderer(rn);
        rw->SetSize(1024, 1024);
        vtkSmartPointer<vtkRenderWindowInteractor> rwi = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        rwi->SetRenderWindow(rw);
        rw->Render();
        rwi->Start();
        
        toCellFilter->SetScalarComponent(0, "alpha.water", 0);
        //toCellFilter->SetScalarComponent(0, "U", 0);
        
        //startRendering(toCellFilter->GetOutputPort());
    } else if (reader->IsFileStructuredGrid()) {
        std::cout << "output is structured grid" << std::endl;
    }
    
    return EXIT_SUCCESS;
}


