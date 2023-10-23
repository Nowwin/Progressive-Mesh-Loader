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


