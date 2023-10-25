#include "MeshMathematics.hpp"

#define BOUNDARY_COST 1.0


bool operator>(const EdgeCollapseTarget &a, const EdgeCollapseTarget &b){
    return a.cost > b.cost;
}

//Helper Methods

//Helper method to find the boundary half edge
HalfEdge* Simplification::FindBoundaryEdgeIncidentToVertexInCW(HalfEdge *baseHalfEdge)
{
    // Try to find a boundary edge by checking incident edges in CW direction
    HalfEdge *hep = baseHalfEdge;
    do {
        // If the current half-edge doesn't have a mate, it's a boundary
        if (hep->mate == nullptr) return hep;

        // Move to the next half-edge in the CW direction
        hep = hep->mate->next;
    } while (hep != baseHalfEdge);

    // Returning the baseHalfEdge if no boundary edge is found.
    // Ideally, this shouldn't happen if the input is guaranteed to have a boundary edge.
    return baseHalfEdge;
}

//updation of quad Error for vertex
void Simplification::CumulateQ(VertexIter &vi, const glm::vec3 &normal, double d) {
    float a = normal.x;
    float b = normal.y;
    float c = normal.z;

    vi->QuadError[0] += a * a;
    vi->QuadError[1] += a * b;
    vi->QuadError[2] += a * c;
    vi->QuadError[3] += a * d;
    vi->QuadError[4] += b * b;
    vi->QuadError[5] += b * c;
    vi->QuadError[6] += b * d;
    vi->QuadError[7] += c * c;
    vi->QuadError[8] += c * d;
    vi->QuadError[9] += d * d;
}

//Creates the Quad error matrxi from the variables
glm::mat4 Simplification::ComputeCombinedQuadric(VertexIter &v0, VertexIter &v1) {
    // Convert the QuadError arrays to glm::mat4
    glm::mat4 q0 = ConvertArrayToMat4(v0->QuadError);
    glm::mat4 q1 = ConvertArrayToMat4(v1->QuadError);

    return q0 + q1; // Return the combined quadric
}

//Just helper to convert to matrix
glm::mat4 Simplification::ConvertArrayToMat4(const float arr[10]) {
    glm::mat4 mat(0.0f);
    mat[0][0] = arr[0];
    mat[0][1] = mat[1][0] = arr[1];
    mat[0][2] = mat[2][0] = arr[2];
    mat[0][3] = mat[3][0] = arr[3];
    mat[1][1] = arr[4];
    mat[1][2] = mat[2][1] = arr[5];
    mat[1][3] = mat[3][1] = arr[6];
    mat[2][2] = arr[7];
    mat[2][3] = mat[3][2] = arr[8];
    mat[3][3] = arr[9];
    return mat;
}

//Functions to solve a linear system
bool Simplification::SolveLinearSystem(const glm::mat4 &matrix, const glm::vec4 &rhs, glm::vec4 &solution) {
    glm::mat4 inverseMatrix = glm::inverse(matrix);
    if (glm::determinant(matrix) == 0.0f) {
        return false; // matrix is singular, cannot solve
    }
    solution = inverseMatrix * rhs;
    return true;
}

//Compute the final cost
float Simplification::ComputeCost(const glm::mat4 &newQ, const glm::vec4 &solution) {
    glm::vec4 temp = newQ * solution;
    return glm::dot(solution, temp);
}

//Assign Quad error for Vertices
void Simplification::AssignInitialQ() {
    for (VertexIter vi = mesh->vertices.begin(); vi != mesh->vertices.end(); vi++)
    {
        //Initialize all quad error amtrix as zero
        for (size_t i = 0; i < 10; i++)
        {
            vi->QuadError[i] = 0.0f;
        }

        HalfEdge* startHalfEdge;
        HalfEdge* endHalfEdge;

        if(!vi->isBoundary) 
            startHalfEdge = vi->neighborHe;
        else                        
            startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(vi->neighborHe);

        HalfEdge* hep = startHalfEdge;

        do {
            //Adding the Quad for each face, vertex, face and the d in n.p+d
            CumulateQ(vi, hep->face->normal_, -glm::dot(hep->face->normal_, hep->face->halfedge[0].vertex->position_));

            if(vi->isBoundary && hep->prev->mate == NULL) {
                endHalfEdge = hep->prev;
                break;
            }

            hep = hep->prev->mate;
        } while(hep != startHalfEdge && hep != NULL);

        //Extra penality for boundaries
        if(vi->isBoundary) {
            // Add pseudo face information to vi->Q

            glm::vec3 boundaryVector = startHalfEdge->next->vertex->position_ - startHalfEdge->vertex->position_;
            glm::vec3 pseudoNormal = glm::cross(boundaryVector, startHalfEdge->face->normal_);
            pseudoNormal = glm::normalize(pseudoNormal);

            CumulateQ(vi, pseudoNormal, -glm::dot(pseudoNormal, startHalfEdge->vertex->position_));

            boundaryVector = endHalfEdge->next->vertex->position_ - endHalfEdge->vertex->position_;
            pseudoNormal = glm::cross(boundaryVector, endHalfEdge->face->normal_);
            pseudoNormal = glm::normalize(pseudoNormal);

            CumulateQ(vi, pseudoNormal, -glm::dot(pseudoNormal, endHalfEdge->vertex->position_));
        }
        
    }
    
}

//Calculates the cost and the cooridnate after collapse for an Edge
void Simplification::ComputeOptimalCoordAndCost(EdgeIter &ei)
{
    VertexIter v0 = ei->halfedge[0]->vertex;
    VertexIter v1 = ei->halfedge[0]->next->vertex;

    glm::mat4 newQ = ComputeCombinedQuadric(v0, v1);
    glm::mat4 matrix;
    glm::vec4 rhs = { 0.0f, 0.0f, 0.0f, 1.0f };
    glm::vec4 solution;

    //Preparing the Quad matrix to be solved
    PrepareMatrix(matrix, newQ);

    float cost;
    glm::vec3 optimalCoord;

    if( SolveLinearSystem(matrix, rhs, solution) )
    {
        cost = ComputeCost(newQ, solution);
        optimalCoord = glm::vec3(solution.x, solution.y, solution.z);
    }
    else
    {
        cost = 0.0f;
        optimalCoord = (v0->isBoundary) ? v0->position_ : v1->position_;
    }

    if(v0->isBoundary || v1->isBoundary) 
        cost += BOUNDARY_COST;

    //Storing the Edge Collapse info in priority Queue based on cost
    heap.push(EdgeCollapseTarget(ei, cost, optimalCoord, ect_id_base));
    ei->ect_id = ect_id_base++;

}

//Checking for residula fins
bool Simplification::IsFinWillNotBeCreated(EdgeIter &edgeIter)
{
    // Retrieve the half-edge associated with the edge we're considering for collapse.
    HalfEdge* collapseEdge = edgeIter->halfedge[0];

    // Retrieve the vertices at the ends of this edge.
    VertexIter startVertex = collapseEdge->vertex;
    VertexIter endVertex = collapseEdge->next->vertex;

    HalfEdge* initialHalfEdgeStart, *initialHalfEdgeEnd;

    // For non-boundary vertices, set the starting half-edge to the collapseEdge.
    // For boundary vertices, find the boundary edge in a clockwise manner.
    initialHalfEdgeStart = (startVertex->isBoundary) ? FindBoundaryEdgeIncidentToVertexInCW(collapseEdge) : collapseEdge;
    initialHalfEdgeEnd = (endVertex->isBoundary) ? FindBoundaryEdgeIncidentToVertexInCW(collapseEdge->next) : collapseEdge->next;

    HalfEdge* currentHalfEdgeStart = initialHalfEdgeStart;

    do {
        HalfEdge* currentHalfEdgeEnd = initialHalfEdgeEnd;

        do {
            // Check if the vertices adjacent to currentHalfEdgeStart and currentHalfEdgeEnd are the same.
            if (currentHalfEdgeStart->next->vertex == currentHalfEdgeEnd->next->vertex ||
                (!currentHalfEdgeEnd->prev->mate && currentHalfEdgeStart->next->vertex == currentHalfEdgeEnd->prev->vertex)) {
                
                VertexIter sharedVertex = currentHalfEdgeStart->next->vertex;

                // If the shared vertex is not one of the two end-points of the collapsing edge, return false.
                if (sharedVertex != collapseEdge->prev->vertex &&
                    (collapseEdge->mate && sharedVertex != collapseEdge->mate->prev->vertex)) {

                    return false;
                }
            }

            // Check for boundary vertices.
            if (!currentHalfEdgeStart->prev->mate) {

                // Check if the previous vertex of currentHalfEdgeStart and the adjacent vertex of currentHalfEdgeEnd are the same.
                if (currentHalfEdgeStart->prev->vertex == currentHalfEdgeEnd->next->vertex ||
                    (!currentHalfEdgeEnd->prev->mate && currentHalfEdgeStart->prev->vertex == currentHalfEdgeEnd->prev->vertex)) {

                    VertexIter sharedBoundaryVertex = currentHalfEdgeStart->prev->vertex;

                    // Again, ensure that if these vertices become neighbors, they won't form a fin.
                    if (sharedBoundaryVertex != collapseEdge->prev->vertex &&
                        (collapseEdge->mate && sharedBoundaryVertex != collapseEdge->mate->prev->vertex)) {

                        return false;
                    }
                }
            }

            // Traverse to the previous connected half-edge.
            currentHalfEdgeEnd = currentHalfEdgeEnd->prev->mate;

        } while (currentHalfEdgeEnd && currentHalfEdgeEnd != initialHalfEdgeEnd);

        // Similarly, traverse to the previous connected half-edge for the start vertex.
        currentHalfEdgeStart = currentHalfEdgeStart->prev->mate;

    } while (currentHalfEdgeStart && currentHalfEdgeStart != initialHalfEdgeStart);

    // If all checks passed, no fin will be created by collapsing this edge.
    return true;
}

//Removes the faces associated with the edge to be collapsed
void Simplification::InactivateFaces(HalfEdge* hepCollapse) {
    //Removing the first face
    hepCollapse->face->isActive = false;

    //Removing the second face if exits
    if (hepCollapse->mate != NULL)
    {
        hepCollapse->mate->face->isActive = false;
        n_active_faces--;
    }

    n_active_faces--;
}

//Stores vertex information for restoration
void Simplification::StoreVertexSplit(EdgeIter &ei, VertexIter &v0, VertexIter &v1) {
    vertexSplitTarget.push(VertexSplitTarget());
    vertexSplitTarget.top().ei = ei;
    vertexSplitTarget.top().v1OriginalCoord = v1->position_;
    vertexSplitTarget.top().v1OriginalIsBoundary = v1->isBoundary;

}

void Simplification::CollectFacesAroundVertices(EdgeIter &ei, VertexIter &v0, VertexIter &v1,  std::vector<FaceIter> &facesOriginallyIncidentToV0OrV1) {
    HalfEdge *hepCollapse = ei->halfedge[0];
    
    //Storing all faces for vertex V0
    HalfEdge *startHalfEdge;
    if(v0->isBoundary == false) {
        startHalfEdge = hepCollapse;
    } else startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(hepCollapse);

    HalfEdge *hep = startHalfEdge;
    do{
        facesOriginallyIncidentToV0OrV1.push_back(hep->face);

        hep = hep->prev->mate;
    }while(hep != startHalfEdge && hep != NULL);

    //Storing all faces for vertex V0
    if(v1->isBoundary == false) startHalfEdge = hepCollapse->next;
    else                        startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(hepCollapse->next);

    hep = startHalfEdge;
    do{
        facesOriginallyIncidentToV0OrV1.push_back(hep->face);

        hep = hep->prev->mate;
    }while(hep != startHalfEdge && hep != NULL);    
}

void Simplification::ReplaceVerticesOfHalfEdges(VertexIter &v0, VertexIter &v1) {
    HalfEdge *startHalfEdge;
    if(!v0->isBoundary) startHalfEdge = v0->neighborHe;
    else                startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(v0->neighborHe);

    HalfEdge *hep = startHalfEdge;
    do {
        if(hep->face->isActive) {
            hep->vertex = v1;
        }

        if(hep->prev->mate == NULL) break;
        hep = hep->prev->mate;
    } while(hep != startHalfEdge && hep != NULL);
}

void Simplification::UpdateEdgeMateInfo(HalfEdge* hepCollapse) {
    hepCollapse->edge->isActive = false;

    if(hepCollapse->next->mate != NULL) hepCollapse->next->mate->mate = hepCollapse->prev->mate;
    if(hepCollapse->prev->mate != NULL) hepCollapse->prev->mate->mate = hepCollapse->next->mate;

    hepCollapse->prev->edge->isActive = false;

    if(hepCollapse->next->edge->halfedge[0] == hepCollapse->next) {
        hepCollapse->next->edge->halfedge[0] = hepCollapse->prev->mate;
    } else {
        hepCollapse->next->edge->halfedge[1] = hepCollapse->prev->mate;
    }

    // Ensuring halfedge[0] is not NULL for edge
    if(hepCollapse->next->edge->halfedge[0] == NULL) {
        if(hepCollapse->next->edge->halfedge[1] != NULL) {
            std::swap(hepCollapse->next->edge->halfedge[0], hepCollapse->next->edge->halfedge[1]);
        } else {
            hepCollapse->next->edge->isActive = false;
        }
    }

    if(hepCollapse->prev->mate != NULL) hepCollapse->prev->mate->edge = hepCollapse->next->edge;

    // Handling the mate side
    if(hepCollapse->mate != NULL) {
        if(hepCollapse->mate->next->mate != NULL) hepCollapse->mate->next->mate->mate = hepCollapse->mate->prev->mate;
        if(hepCollapse->mate->prev->mate != NULL) hepCollapse->mate->prev->mate->mate = hepCollapse->mate->next->mate;

        hepCollapse->mate->next->edge->isActive = false;

        if(hepCollapse->mate->prev->edge->halfedge[0] == hepCollapse->mate->prev) {
            hepCollapse->mate->prev->edge->halfedge[0] = hepCollapse->mate->next->mate;
        } else {
            hepCollapse->mate->prev->edge->halfedge[1] = hepCollapse->mate->next->mate;
        }

        // Ensuring halfedge[0] is not NULL for the mate's edge
        if(hepCollapse->mate->prev->edge->halfedge[0] == NULL) {
            if(hepCollapse->mate->prev->edge->halfedge[1] != NULL) {
                std::swap(hepCollapse->mate->prev->edge->halfedge[0], hepCollapse->mate->prev->edge->halfedge[1]);
            } else {
                hepCollapse->mate->prev->edge->isActive = false;
            }
        }

        if(hepCollapse->mate->next->mate != NULL) hepCollapse->mate->next->mate->edge = hepCollapse->mate->prev->edge;
    }
}

void Simplification::UpdateVertexNormal(VertexIter &v) {
    glm::vec3 normalSum(0.0f, 0.0f, 0.0f);
    int faceCount = 0;

    HalfEdge *startHalfEdge;
    if(!v->isBoundary) startHalfEdge = v->neighborHe;
    else                startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(v->neighborHe);

    HalfEdge *hep = startHalfEdge;
    do {
        if(hep->face->isActive) {
            normalSum += hep->face->normal_;  // Assuming you have face normals precomputed or available
            faceCount++;
        }

        if(hep->prev->mate == NULL) break;
        hep = hep->prev->mate;
    } while(hep != startHalfEdge && hep != NULL);

    if(faceCount > 0) {
        v->normal_ = glm::normalize(normalSum / static_cast<float>(faceCount));
    }
}


void Simplification::RemoveEdge(EdgeIter &ei, glm::vec3 optimalCoord, bool isFirstCollapse) {
    HalfEdge *hepCollapse = ei->halfedge[0];
    VertexIter v0 = hepCollapse->vertex;
    VertexIter v1 = hepCollapse->next->vertex;
    InactivateFaces(hepCollapse);
    StoreVertexSplit(ei, v0, v1);

    std::vector<FaceIter> facesOriginallyIncidentToV0OrV1;
    CollectFacesAroundVertices(ei, v0, v1, facesOriginallyIncidentToV0OrV1);
    ReplaceVerticesOfHalfEdges(v0, v1);
    UpdateEdgeMateInfo(hepCollapse);   
    UpdateVertexNormal(v1); 
}


//Process edges based on cost from the heap
bool Simplification::ProcessEdgeCollapseHeap() {
    while(!heap.empty()) {
        EdgeCollapseTarget ect = heap.top();
        heap.pop();

        if(ect.ei->isActive == true && ect.id == ect.ei->ect_id) {
            if(IsFinWillNotBeCreated(ect.ei)) {
                RemoveEdge(ect.ei, ect.optimalCoord, true); 
                return true;
            }
            else {
                suspendedEdgeCollapseTarget.push_back(ect);
            }
        }
    }
    return false;
}

bool Simplification::ProcessReaddedEdgeCollapseTarget() {
    if (readdedEdgeCollapseTarget.empty() == false) {
        RemoveEdge(readdedEdgeCollapseTarget.top().ei, readdedEdgeCollapseTarget.top().optimalCoord, false);
        readdedEdgeCollapseTarget.pop();
        return true;
    }
    return false;
}

bool Simplification::ProcessSuspendedEdgeCollapseTarget() {
    std::list<EdgeCollapseTarget>::iterator ecti = suspendedEdgeCollapseTarget.begin();
    while (ecti != suspendedEdgeCollapseTarget.end()) {
        
        if(ecti->ei->isActive == false || ecti->id != ecti->ei->ect_id){
            // obsolete. delete this
            ecti = suspendedEdgeCollapseTarget.erase(ecti);
        } else {
            if( IsFinWillNotBeCreated(ecti->ei) ){
                RemoveEdge(ecti->ei, ecti->optimalCoord, true);
                ecti = suspendedEdgeCollapseTarget.erase(ecti);
                return true;
            } else {
                ecti++;
            }
        }
    }
    return false;
}




//Actual Functions ----------------------

//------------Initialization-------------------------------------

//Initialization
void Simplification::InitSimplification(Mesh *mesh_in) {
    this->mesh = mesh_in;
    this->n_active_faces = this->mesh->n_faces;

    AssignInitialQ();

    for(EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei)
        ComputeOptimalCoordAndCost(ei);
}

//-----------------------------------------------------------------


//------------Mesh Simplification: Collapsing an Edge-------------------------------------

bool Simplification::EdgeCollapse() {
    if (n_active_faces < 3) return false;
    if (ProcessReaddedEdgeCollapseTarget()) return true;
    if (ProcessSuspendedEdgeCollapseTarget()) return true;
    return ProcessEdgeCollapseHeap();
}

//-----------------------------------------------------------------