//
//  Tetragrid.cpp
//  gridReader
//
//  Created by Franz Neubert on 05/02/16.
//  Copyright © 2016 Franz Neubert. All rights reserved.
//

#include <set>
#include <vector>
#include "Calculator.h"
#include "Helper.hpp"
#include "Tetragrid.hpp"
#include "Edge.hpp"

using namespace std;

Tetragrid *Tetragrid::createGrid(vtkUnstructuredGrid *grid) {
    vector<Cell*> cells;
    vector<Vertex*> vertices = Vertex::verticesFromGrid(grid);
    vector<Edge*> edges;
    
    
    for (vtkIdType cell = 0; cell < grid->GetNumberOfCells(); cell++) {
        vector<Edge*> cellEdges;
        for (vtkIdType edge = 0; edge < grid->GetCell(cell)->GetNumberOfEdges(); edge++) {
            int vtkA = (int) grid->GetCell((int)cell)->GetEdge((int)edge)->GetPointId(0);
            int vtkB = (int) grid->GetCell((int)cell)->GetEdge((int)edge)->GetPointId(1);
            Edge *checkEdge = Edge::isEdge(vertices.at(vtkA), vertices.at(vtkB), &edges);
            if (checkEdge == nullptr) {
                checkEdge = new Edge(vertices.at(vtkA), vertices.at(vtkB));
                edges.push_back(checkEdge);
            }
            cellEdges.push_back(checkEdge);
        }
        vector<Vertex*> cellVertices;
        for (vtkIdType point = 0; point < grid->GetCell(cell)->GetNumberOfPoints(); point++) {
            cellVertices.push_back(vertices.at(grid->GetCell(cell)->GetPointId((int)point)));
        }
        cells.push_back(new Cell((int)cell, cellVertices, cellEdges));
    }
    
    Tetragrid *tetragrid = new Tetragrid(cells, edges, vertices);
    return tetragrid;
}

void Tetragrid::precalculations() {
    //this->precalcCells();
    this->precalcEdges();
    this->precalceVertices();
}

void Tetragrid::precalcEdges() {
    for (auto edge : this->edges) {
        double coordsA[3];
        double coordsB[3];
        edge->getA()->getCoords(coordsA);
        edge->getB()->getCoords(coordsB);
        edge->setLength(Calculator::calcEdgeLength(coordsA, coordsB));
    }
}

void Tetragrid::precalceVertices() {
    for (auto cell : this->cells) {
        for (auto vertex : cell->vertices) {
            vertex->incidents.push_back(cell);
        }
    }
}