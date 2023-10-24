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

//Initialization
void Simplification::InitSimplification(Mesh *mesh_in) {
    this->mesh = mesh_in;
    this->n_active_faces = this->mesh->n_faces;

    AssignInitialQ();

    for(EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei)
        ComputeOptimalCoordAndCost(ei);
}

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