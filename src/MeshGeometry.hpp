#pragma once

#include<iostream>
#include<vector>
#include <list>
#include<unordered_map>

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

    float QuadError[10]; //This represents the Quad error, all the constants of matrix

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
private:
    std::unordered_map<int, int> vertexIndexMap;    
protected:
    //Reads from an OBJ file - ignores texture and normals
    bool ReadOBJFile(char *filename);
    //Setups half edge and associated face, edge, vertex
    void AddEdgeInfo();
    //Setups the halfedge data for a face
    void MakeCircularList(FaceIter &fi);

public:
    
    std::list<Vertex> vertices;
    std::list<Face>   faces;
    std::list<Edge>   edges;

    int n_vertices, n_faces, n_edges;


    Mesh(){
        n_vertices = n_faces = n_edges = 0;
    }

    //Begins the construction of Mesh
    bool ConstructMeshDataStructure(char *filename);
    //Calculates Face Normal
    void AssignFaceNormal(FaceIter &fi);
    //Calulcate Vertex Normal
    void AssignVertexNormal(VertexIter &vi);
    //Gets the vertex data
    std::vector<GLfloat> GetVertexData();
    //Get the EBO
    std::vector<GLuint> GetIndexData();

};
