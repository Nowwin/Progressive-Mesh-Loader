#include "MeshMathematics.hpp"

#define BOUNDARY_COST 1.0


bool operator>(const EdgeCollapseTarget &a, const EdgeCollapseTarget &b){
    return a.cost > b.cost;
}

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

void Simplification::FindNeighborHalfEdge(VertexIter &v1, std::vector<FaceIter> &facesOriginallyIncidentToV0OrV1)
{
    // Loop through the list of faces originally incident to either v0 or v1
    for(auto &face : facesOriginallyIncidentToV0OrV1)
    {
        // If the face is still active
        if(face->isActive)
        {
            // Loop through the half-edges of the face
            for(int j = 0; j < 3; j++)
            {
                // If the vertex associated with the current half-edge is v1
                if(face->halfedge[j].vertex == v1) 
                {
                    // Set the neighbor half-edge of v1
                    v1->neighborHe = &(face->halfedge[j]);
                    return; // Exit once found
                }
            }
        }
    }
}


//Calculation of quad Error
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

void Simplification::AssignInitialQ() {
    for(VertexIter vi = mesh->vertices.begin(); vi != mesh->vertices.end(); vi++) {

        for(int i = 0; i < 10; i++) vi->QuadError[i] = 0.0f;

        HalfEdge* startHalfEdge;
        HalfEdge* endHalfEdge;

        if(!vi->isBoundary) 
            startHalfEdge = vi->neighborHe;
        else                        
            startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(vi->neighborHe);

        HalfEdge* hep = startHalfEdge;
        do {
            CumulateQ(vi, hep->face->normal_, -glm::dot(hep->face->normal_, hep->face->halfedge[0].vertex->position_));

            if(vi->isBoundary && hep->prev->mate == NULL) {
                endHalfEdge = hep->prev;
                break;
            }

            hep = hep->prev->mate;
        } while(hep != startHalfEdge && hep != NULL);

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

void Simplification::InitSimplification(Mesh *mesh_in)
{
    mesh = mesh_in;

    n_active_faces = mesh->n_faces;

    AssignInitialQ();

    for(EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei)
        ComputeOptimalCoordAndCost(ei);
}

void Simplification::ComputeOptimalCoordAndCost(EdgeIter &ei)
{
    VertexIter v0 = ei->halfedge[0]->vertex;
    VertexIter v1 = ei->halfedge[0]->next->vertex;

    glm::mat4 newQ = ComputeCombinedQuadric(v0, v1);
    glm::mat4 matrix;
    glm::vec4 rhs = { 0.0f, 0.0f, 0.0f, 1.0f };
    glm::vec4 solution;

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

    heap.push(EdgeCollapseTarget(ei, cost, optimalCoord, ect_id_base));
    ei->ect_id = ect_id_base++;

}

bool Simplification::SolveLinearSystem(const glm::mat4 &matrix, const glm::vec4 &rhs, glm::vec4 &solution) {
    glm::mat4 inverseMatrix = glm::inverse(matrix);
    if (glm::determinant(matrix) == 0.0f) {
        return false; // matrix is singular, cannot solve
    }
    solution = inverseMatrix * rhs;
    return true;
}

glm::mat4 Simplification::ComputeCombinedQuadric(VertexIter &v0, VertexIter &v1) {
    // Convert the QuadError arrays to glm::mat4
    glm::mat4 q0 = ConvertArrayToMat4(v0->QuadError);
    glm::mat4 q1 = ConvertArrayToMat4(v1->QuadError);

    return q0 + q1; // Return the combined quadric
}


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

void Simplification::PrepareMatrix(glm::mat4 &matrix, const glm::mat4 &newQ) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            matrix[i][j] = newQ[i][j];
    matrix[3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
}


float Simplification::ComputeCost(const glm::mat4 &newQ, const glm::vec4 &solution) {
    glm::vec4 temp = newQ * solution;
    return glm::dot(solution, temp);
}


