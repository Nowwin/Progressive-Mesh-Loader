#pragma once

#include <queue>
#include <stack>
#include <list>

#include "glad/glad.h"
#include "glm/glm.hpp"

#include "MeshGeometry.hpp"

struct EdgeCollapseTarget {
    EdgeIter ei;
    float cost;
    glm::vec3 optimalCoord;
    int id;

    EdgeCollapseTarget(){}
    EdgeCollapseTarget(EdgeIter &ei_in, float cost_in, const glm::vec3 &optimalCoord_in, int id_in) 
        : ei(ei_in), cost(cost_in), optimalCoord(optimalCoord_in), id(id_in) {}
};


struct VertexSplitTarget {
    EdgeIter ei;
    glm::vec3 v1OriginalCoord;
    bool     v1OriginalIsBoundary;
    std::vector<HalfEdge*> halfedgesAroundV0;

    int id;

    VertexSplitTarget(){}
};

class Simplification {
    Mesh *mesh;

    std::priority_queue <EdgeCollapseTarget, std::deque<EdgeCollapseTarget>, std::greater<EdgeCollapseTarget>> heap;
    std::list<EdgeCollapseTarget> suspendedEdgeCollapseTarget;

    std::stack<VertexSplitTarget>  vertexSplitTarget;
    std::stack<EdgeCollapseTarget> readdedEdgeCollapseTarget;

    int ect_id_base;
    int n_active_faces;

    void AssignInitialQ();
    void CumulateQ(VertexIter &vi, const glm::vec3 &normal, double d);
    void ComputeOptimalCoordAndCost(EdgeIter &ei);
    HalfEdge* FindBoundaryEdgeIncidentToVertexInCW(HalfEdge *baseHalfEdge);
    void FindNeighborHalfEdge(VertexIter &v1, std::vector<FaceIter> &facesOriginallyIncidentToV0OrV1);
    bool IsFinWillNotBeCreated(EdgeIter &ei);
    void RemoveEdge(EdgeIter &ei, glm::vec3 optimalCoord, bool isFirstCollapse);
    glm::mat4 ComputeCombinedQuadric(VertexIter &v0, VertexIter &v1);
    void PrepareMatrix(glm::mat4 &matrix, const glm::mat4 &newQ);
    float ComputeCost(const glm::mat4 &newQ, const glm::vec4 &solution);
    bool SolveLinearSystem(const glm::mat4 &matrix, glm::vec4 &rhs, glm::vec4 &solution);
    glm::mat4 ConvertArrayToMat4(const float arr[10]); 
    bool ProcessReaddedEdgeCollapseTarget();
    bool ProcessSuspendedEdgeCollapseTarget();
    bool ProcessEdgeCollapseHeap();
    void InactivateFaces(HalfEdge* hepCollapse);
    void StoreVertexSplit(EdgeIter &ei, VertexIter &v0, VertexIter &v1);
    void ReplaceVerticesOfHalfEdges(VertexIter &v0, VertexIter &v1);
    void UpdateEdgeMateInfo(HalfEdge* hepCollapse);
    void CollectFacesAroundVertices(EdgeIter &ei, VertexIter &v0, VertexIter &v1,  std::vector<FaceIter> &facesOriginallyIncidentToV0OrV1);
    void UpdateNormalsAroundVertex(VertexIter &v_target);
public:
    Simplification(){ ect_id_base = 0; }

    void InitSimplification(Mesh *mesh_in);
    bool EdgeCollapse();
    void VertexSplit();
    void ControlLevelOfDetail(int step);
    Mesh* GetModifiedMesh() const {
        return mesh;
    }
};


bool validateMesh(Mesh* mesh);