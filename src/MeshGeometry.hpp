#pragma once

#include <list>
#include "glad/glad.h"
#include "glm/glm.hpp"

//This is the interface for Half Edge data structure


struct HalfEdge;
struct Vertex;
struct Face;
struct Edge;

typedef std::list<Vertex>::iterator  VertexIter;
typedef std::list<Face>::iterator    FaceIter;
typedef std::list<Edge>::iterator    EdgeIter;


struct HalfEdge {
    HalfEdge *next, *prev;
    HalfEdge *mate;

    VertexIter vertex;
    FaceIter   face;
    EdgeIter   edge;

    HalfEdge(){
        next = prev = NULL;
        mate = NULL;
    }
};

struct Edge {
    HalfEdge *halfedge[2];
    int id;

    int ect_id; 

    bool isActive;

    Edge(){}
    Edge(HalfEdge *_halfedgeForward, HalfEdge *halfedgeBackward, int n){
        halfedge[0] = _halfedgeForward; 
        halfedge[1] = halfedgeBackward;
        id = n;
        isActive = true;
    };
};


struct Vertex {
    glm::vec3 position_;
    glm::vec3 normal_;
    int id;
    
    HalfEdge *neighborHe;

    bool isBoundary;
    bool isActive;

    double QuadError[10]; //This represents the Quad error, all the constants of matrix

    Vertex(){};
    Vertex(const glm::vec3& position, int n) : position_(position), id(n) {
        neighborHe = NULL;
        isBoundary = false;
        isActive = true;
    }
};


struct Face {
    
    HalfEdge   halfedge[3]; 

    glm::vec3 normal_;
    double area;
    int id;

    bool isActive;

    Face(){}
    Face(VertexIter v0, VertexIter v1, VertexIter v2, int n){
        halfedge[0].vertex = v0;
        halfedge[1].vertex = v1;
        halfedge[2].vertex = v2;
        id = n;
        isActive = true;
    }
};


class Mesh {
protected:

    bool ReadFile(char *filename);
    void AddEdgeInfo();
    void MakeCircularList(FaceIter &fi);

public:
    
    std::list<Vertex> vertices;
    std::list<Face>   faces;
    std::list<Edge>   edges;

    int n_vertices, n_faces, n_edges;


    Mesh(){
        n_vertices = n_faces = n_edges = 0;
    }

    bool ConstructMeshDataStructure(char *filename);
    void AssignFaceNormal(FaceIter &fi);
    void AssignVertexNormal(VertexIter &vi);
    void Display(int mode);

};