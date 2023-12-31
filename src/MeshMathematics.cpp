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
    return hep;
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
bool Simplification::SolveLinearSystem(const glm::mat4 &matrix, glm::vec4 &rhs, glm::vec4 &solution) {
    glm::mat4 mat = matrix;  // Create a copy for manipulation

    // Gaussian elimination
    for (int i = 0; i < 4; i++) {
        // Find pivot
        double maxAbsVal = -1.0;
        int pivotIndex = i;
        for (int j = i; j < 4; j++) {
            if (fabs(mat[j][i]) > maxAbsVal) {
                maxAbsVal = fabs(mat[j][i]);
                pivotIndex = j;
            }
        }

        // Check for near-singularity
        if (maxAbsVal < 1e-6) {
            return false;
        }

        // Swap rows
        std::swap(mat[i], mat[pivotIndex]);
        std::swap(rhs[i], rhs[pivotIndex]);

        // Eliminate
        for (int j = i + 1; j < 4; j++) {
            double factor = mat[j][i] / mat[i][i];
            for (int k = i; k < 4; k++) {
                mat[j][k] -= factor * mat[i][k];
            }
            rhs[j] -= factor * rhs[i];
        }
    }

    // Back substitution
    for (int i = 3; i >= 0; i--) {
        solution[i] = rhs[i];
        for (int j = i + 1; j < 4; j++) {
            solution[i] -= mat[i][j] * solution[j];
        }
        solution[i] /= mat[i][i];
    }

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

void Simplification::PrepareMatrix(glm::mat4 &matrix, const glm::mat4 &newQ) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            matrix[i][j] = newQ[i][j];
    matrix[3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
}

//For debugging
void printMat4(const glm::mat4 &mat) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cerr << mat[i][j] << " ";
        }
        std::cerr << std::endl;
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

    //std::cerr << "My matrix:\n";
    //printMat4(matrix);

    float cost;
    glm::vec3 optimalCoord;

    if( SolveLinearSystem(matrix, rhs, solution) )
    {
        //std::cerr << "The solution is: " << solution.x << " " << solution.y << " " << solution.z << " " << std::endl;
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
            vertexSplitTarget.top().halfedgesAroundV0.push_back(hep);
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

void Simplification::FindNeighborHalfEdge(VertexIter &v1, std::vector<FaceIter> &facesOriginallyIncidentToV0OrV1)
{
    // choose neighborHe of v1 from a face that is still active

    for(unsigned int i = 0; i < facesOriginallyIncidentToV0OrV1.size(); i++){
        if(facesOriginallyIncidentToV0OrV1[i]->isActive == true){
            for(int j = 0; j < 3; j++){
                if(facesOriginallyIncidentToV0OrV1[i]->halfedge[j].vertex == v1) {
                    v1->neighborHe = &(facesOriginallyIncidentToV0OrV1[i]->halfedge[j]);
                    break;
                }
            }
            break;
        }
    }
}

/*void Simplification::RemoveEdge(EdgeIter &ei, glm::vec3 optimalCoord, bool isFirstCollapse) {
    HalfEdge *hepCollapse = ei->halfedge[0];
    VertexIter v0 = hepCollapse->vertex;
    VertexIter v1 = hepCollapse->next->vertex;
    InactivateFaces(hepCollapse);
    StoreVertexSplit(ei, v0, v1);

    std::vector<FaceIter> facesOriginallyIncidentToV0OrV1;
    CollectFacesAroundVertices(ei, v0, v1, facesOriginallyIncidentToV0OrV1);
    ReplaceVerticesOfHalfEdges(v0, v1);
    FindNeighborHalfEdge(v1, facesOriginallyIncidentToV0OrV1);

    v1->position_ = optimalCoord;
    if(isFirstCollapse){
        // add v0's "Q" to v1's "Q"
        for(int i = 0; i < 10; i++) v1->QuadError[i] += v0->QuadError[i];
    }

    UpdateEdgeMateInfo(hepCollapse);   
    UpdateVertexNormal(v1); 
}*/

void Simplification::RemoveEdge(EdgeIter &ei, glm::vec3 optimalCoord, bool isFirstCollapse)
{
    HalfEdge *hepCollapse = ei->halfedge[0];
    //std::cerr << ei->halfedge[0]->edge->id << std::endl;
    //std::cerr << ei->halfedge[1]->edge->id << std::endl;
    HalfEdge *startHalfEdge = nullptr;
    HalfEdge * hep = startHalfEdge;

    VertexIter v0 = hepCollapse->vertex;
    VertexIter v1 = hepCollapse->next->vertex;

    //std::cout << "Positions to collapse:\n" << std::endl;
    //std::cout << v0->position_.x << ", " << v0->position_.y << ", " << v0->position_.z << std::endl;
    //std::cout << v1->position_.x << ", " << v1->position_.y << ", " << v1->position_.z << std::endl; 
   
    InactivateFaces(hepCollapse);
    StoreVertexSplit(ei, v0, v1);
    

    std::vector<FaceIter> facesOriginallyIncidentToV0OrV1; 
    CollectFacesAroundVertices(ei, v0, v1, facesOriginallyIncidentToV0OrV1);
    ReplaceVerticesOfHalfEdges(v0, v1);

    //Aligning v1 to the new coordinate and updating quad
    v1->position_ = optimalCoord;

    if(isFirstCollapse){
        for(int i = 0; i < 10; i++) v1->QuadError[i] += v0->QuadError[i];
    }

    //Half Edge and Edge Updates
    hepCollapse->edge->isActive = false;

    if(hepCollapse->next->mate != NULL) hepCollapse->next->mate->mate = hepCollapse->prev->mate;
    if(hepCollapse->prev->mate != NULL) hepCollapse->prev->mate->mate = hepCollapse->next->mate;

    hepCollapse->prev->edge->isActive = false;

    if(hepCollapse->next->edge->halfedge[0] == hepCollapse->next) 
        hepCollapse->next->edge->halfedge[0] = hepCollapse->prev->mate;
    else                                         
        hepCollapse->next->edge->halfedge[1] = hepCollapse->prev->mate;

    // edge->halfedge[0] should not be NULL
    if(hepCollapse->next->edge->halfedge[0] == NULL){
        if(hepCollapse->next->edge->halfedge[1] != NULL){
            // swap
            hepCollapse->next->edge->halfedge[0] = hepCollapse->next->edge->halfedge[1];
            hepCollapse->next->edge->halfedge[1] = NULL;
        }else{ // hepCollapse->next->edge->halfedge[0] == NULL and hepCollapse->next->edge->halfedge[1] == NULL
            // edge becomes degenerate
            hepCollapse->next->edge->isActive = false;
        }
    }

    if(hepCollapse->prev->mate != NULL) hepCollapse->prev->mate->edge = hepCollapse->next->edge;


    if(hepCollapse->mate != NULL){ // when "ect->ei" is not a boundary edge
        if(hepCollapse->mate->next->mate != NULL) hepCollapse->mate->next->mate->mate = hepCollapse->mate->prev->mate;
        if(hepCollapse->mate->prev->mate != NULL) hepCollapse->mate->prev->mate->mate = hepCollapse->mate->next->mate;

        hepCollapse->mate->next->edge->isActive = false;

        if(hepCollapse->mate->prev->edge->halfedge[0] == hepCollapse->mate->prev) 
            hepCollapse->mate->prev->edge->halfedge[0] = hepCollapse->mate->next->mate;
        else                                                      
            hepCollapse->mate->prev->edge->halfedge[1] = hepCollapse->mate->next->mate;

        // edge->halfedge[0] should not be NULL
        if(hepCollapse->mate->prev->edge->halfedge[0] == NULL){
            if(hepCollapse->mate->prev->edge->halfedge[1] != NULL){
                // swap
                hepCollapse->mate->prev->edge->halfedge[0] = hepCollapse->mate->prev->edge->halfedge[1];
                hepCollapse->mate->prev->edge->halfedge[1] = NULL;
            }else{ // hepCollapse->mate->prev->edge->halfedge[0] == NULL and hepCollapse->mate->prev->edge->halfedge[1] == NULL
                // edge becomes degenerate
                hepCollapse->mate->prev->edge->isActive = false;
            }
        }

        if(hepCollapse->mate->next->mate != NULL) hepCollapse->mate->next->mate->edge = hepCollapse->mate->prev->edge;
    
    }

   if( hepCollapse->next->edge->isActive == false && 
       (hepCollapse->mate == NULL || hepCollapse->mate->prev->edge->isActive == false) ){
       // face disappeared after edge collapsing
       return;
   }



    if(v0->isBoundary) v1->isBoundary = true;

    FindNeighborHalfEdge(v1, facesOriginallyIncidentToV0OrV1);

    if(v1->isBoundary == false) startHalfEdge = v1->neighborHe;
    else                        startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(v1->neighborHe);

    // update cost and optimal vertex coordinate of incident edges, and normals of incident faces, and incident vertices' neighborHe;
    hep = startHalfEdge;
    do{

       if(isFirstCollapse) ComputeOptimalCoordAndCost(hep->edge);

        mesh->AssignFaceNormal(hep->face);

        hep->next->vertex->neighborHe = hep->next;

        if(hep->prev->mate == NULL){
            if(isFirstCollapse) ComputeOptimalCoordAndCost(hep->prev->edge);    
            hep->prev->vertex->neighborHe = hep->prev;
            break;   
        }

        hep = hep->prev->mate;
    }while(hep != startHalfEdge && hep != NULL);


    // Updating the Normals
    mesh->AssignVertexNormal(v1);

    //std::cerr << ei->halfedge[0]->edge->id << std::endl;
    //std::cerr << ei->halfedge[1]->edge->id << std::endl;

    hep = startHalfEdge;
    do{
        mesh->AssignVertexNormal(hep->next->vertex);

        if(hep->prev->mate == NULL){ 
            mesh->AssignVertexNormal(hep->prev->vertex);
            break;   
        }

        hep = hep->prev->mate;
    }while(hep != startHalfEdge);

    //std::cerr << ei->halfedge[0]->edge->id << std::endl;
    //std::cerr << ei->halfedge[1]->edge->id << std::endl;
}



//Process edges based on cost from the heap
bool Simplification::ProcessEdgeCollapseHeap() {
    while(!heap.empty()) {
        EdgeCollapseTarget ect = heap.top();
        heap.pop();
        //std::cerr << "loop-check\n";
        if(ect.ei->isActive == true && ect.id == ect.ei->ect_id) {
            if(IsFinWillNotBeCreated(ect.ei)) {
                //std::cerr << "loop-check-1\n";
                //std::cerr << ect.optimalCoord.x << " " << ect.optimalCoord.y << " " << ect.optimalCoord.z << std::endl;
                RemoveEdge(ect.ei, ect.optimalCoord, true);
                //std::cerr << "loop-check-2\n"; 
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


void Simplification::UpdateNormalsAroundVertex(VertexIter &v_target) {
    // Determine the starting half-edge
    HalfEdge *startHalfEdge;
    if(v_target->isBoundary == false) {
        startHalfEdge = v_target->neighborHe;
    } else {
        startHalfEdge = FindBoundaryEdgeIncidentToVertexInCW(v_target->neighborHe);
    }

    // Update normals of incident faces
    HalfEdge *hep = startHalfEdge;
    do {
        mesh->AssignFaceNormal(hep->face);
        hep = hep->prev->mate;
    } while (hep != startHalfEdge && hep != NULL);

    // Update the normal of the target vertex
    mesh->AssignVertexNormal(v_target);

    // Update normals of adjacent vertices
    hep = startHalfEdge;
    do {
        mesh->AssignVertexNormal(hep->next->vertex);
        if (hep->prev->mate == NULL) {
            mesh->AssignVertexNormal(hep->prev->vertex);
            break;
        }
        hep = hep->prev->mate;
    } while (hep != startHalfEdge);
}


//Error checking

bool isValidLoop(HalfEdge* startHalfEdge) {
    HalfEdge* currentHalfEdge = startHalfEdge;
    const int MAX_COUNT = 10;  // assuming triangular faces for simplicity
    int count = 0;

    do {
        if (currentHalfEdge->next->prev != currentHalfEdge) return false;
        if (currentHalfEdge->prev->next != currentHalfEdge) return false;

        currentHalfEdge = currentHalfEdge->next;
        count++;
    } while (currentHalfEdge != startHalfEdge && count < MAX_COUNT);

    return currentHalfEdge == startHalfEdge;
}

bool validateMesh(Mesh* mesh) {
    

    // Validate half-edges
    for (EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei) {
        if (ei->halfedge[0] == nullptr || !isValidLoop(ei->halfedge[0])) {
            std::cerr << "Half Edge issue" << std::endl;
            return false;
        }
        if (ei->halfedge[1] != nullptr && !isValidLoop(ei->halfedge[1])) {
            std::cerr << "Half Edge Mate issue" << std::endl;
            return false;
        }
    }

    // Validate edges
    for (EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei) {
        if (ei->isActive && (&(*(ei->halfedge[0]->edge)) != &(*ei))) {
            std::cerr << ei->halfedge[0]->edge->id << std::endl;
            std::cerr << ei->halfedge[1]->edge->id << std::endl;
            std::cerr << ei->id << std::endl;

            std::cerr << "Half Edge in Edge issue" << std::endl;
            return false;
        } 
        if (ei->isActive && ei->halfedge[1] != nullptr && (&(*ei->halfedge[1]->edge)) != &(*ei)) {
            std::cerr << ei->halfedge[0]->edge->id << std::endl;
            std::cerr << ei->halfedge[1]->edge->id << std::endl;
            std::cerr << ei->id << std::endl;
            std::cerr << "Half Edge Mate in Edge issue" << std::endl;
            return false;
        }
    }

    // Validate faces
    for (FaceIter fi = mesh->faces.begin(); fi != mesh->faces.end(); ++fi) {
        for (int i = 0; i < 3; ++i) {
            if (!isValidLoop(&fi->halfedge[i])) {
                std::cerr << "Face issue" << std::endl;
                return false;
            } 
        }
    }

    return true;
}



//Actual Functions ----------------------

//------------Initialization-------------------------------------

//Initialization
void Simplification::InitSimplification(Mesh *mesh_in) {
    std::cerr << "Intialization started\n";
    this->mesh = mesh_in;
    this->n_active_faces = this->mesh->n_faces;

    AssignInitialQ();
    std::cerr << "Intialization of Q done\n";

    for(EdgeIter ei = mesh->edges.begin(); ei != mesh->edges.end(); ++ei)
        ComputeOptimalCoordAndCost(ei);
    
    std::cerr << "Intialization Finished\n";
}

//-----------------------------------------------------------------


//------------Mesh Simplification: Collapsing an Edge-------------------------------------

bool Simplification::EdgeCollapse() {
    if (n_active_faces < 3) return false;
    //std::cerr << "check - 1\n";
    if (ProcessReaddedEdgeCollapseTarget()) return true;
    //std::cerr << "check - 2\n";
    if (ProcessSuspendedEdgeCollapseTarget()) return true;
    //std::cerr << "check - 3\n";
    return ProcessEdgeCollapseHeap();
}

//-----------------------------------------------------------------

//------------Mesh Simplification: Restoring an Edge-------------------------------------

void Simplification::VertexSplit() {
    if(vertexSplitTarget.empty() == false){

        //Identifying one of the half-edge from the edge that was collapsed
        HalfEdge *hepCollapsed = vertexSplitTarget.top().ei->halfedge[0];

        //std::cerr << vertexSplitTarget.top().ei->halfedge[0]->edge->id << std::endl;
        //std::cerr << vertexSplitTarget.top().ei->halfedge[1]->edge->id << std::endl;
        //std::cerr << hepCollapsed->edge->id << std::endl;
        //std::cerr << hepCollapsed->prev->edge->id << std::endl;
        //std::cerr << hepCollapsed->prev->mate->edge->id << std::endl;
        //Start and end point of the half edge
        VertexIter v0 = hepCollapsed->vertex; 
        VertexIter v1 = hepCollapsed->next->vertex; 

        //std::cout << "Positions to restore:\n" << std::endl;
        //std::cout << v0->position_.x << ", " << v0->position_.y << ", " << v0->position_.z << std::endl;
        //std::cout << v1->position_.x << ", " << v1->position_.y << ", " << v1->position_.z << std::endl;
        //std::cout << "Orginal Position to restore:\n" << std::endl; 
        //std::cout << vertexSplitTarget.top().v1OriginalCoord.x << ", " << vertexSplitTarget.top().v1OriginalCoord.y
        // << ", " << vertexSplitTarget.top().v1OriginalCoord.z << std::endl;

        // re-add edge collapse, no need to calculate anything since they were removed in order
        readdedEdgeCollapseTarget.push( EdgeCollapseTarget(vertexSplitTarget.top().ei, -1.0f, v1->position_, -1) );

        //Restoring the original coordinate
        v1->position_ = vertexSplitTarget.top().v1OriginalCoord;
        v1->isBoundary = vertexSplitTarget.top().v1OriginalIsBoundary;

        //Edge restoration
        //We only need to restore the previous and collapsed edge explicitly
        //v1 was just just shifted in corrdinates but still have v1 to vx edge
        hepCollapsed->edge->isActive = true;

        FaceIter currentFace = hepCollapsed->face;
        currentFace->isActive = true;
        n_active_faces++;

        //Ensuring the theoretical consistancy for mate and mate's
        if(hepCollapsed->next->mate != NULL) {
            hepCollapsed->next->mate->mate = hepCollapsed->next;
            hepCollapsed->next->mate->edge = hepCollapsed->next->edge;
        } 
        if(hepCollapsed->prev->mate != NULL) {
            hepCollapsed->prev->mate->mate = hepCollapsed->prev;
            hepCollapsed->prev->mate->edge = hepCollapsed->prev->edge;
        } 

        //Restoring the previous edge
        hepCollapsed->prev->edge->isActive = true;
        

        if(hepCollapsed->next->edge->halfedge[0] == hepCollapsed->next->mate) 
            hepCollapsed->next->edge->halfedge[1] = hepCollapsed->next;
        else {
            hepCollapsed->next->edge->halfedge[0] = hepCollapsed->next;
            /*if (hepCollapsed->next->edge->halfedge[1] != nullptr)
            {
                hepCollapsed->next->edge->halfedge[1] = hepCollapsed->next->mate;
            }*/
            
        } 

        if(hepCollapsed->prev->edge->halfedge[0] == hepCollapsed->prev->mate) 
            hepCollapsed->prev->edge->halfedge[1] = hepCollapsed->prev;
        else {
            hepCollapsed->prev->edge->halfedge[0] = hepCollapsed->prev;
            if (hepCollapsed->prev->edge->halfedge[1] != nullptr)
            {
                hepCollapsed->prev->edge->halfedge[1] = hepCollapsed->prev->mate;
            }
            
        }                 
            
  
        //Restoring Face's vertex half edge connection
        for(int i = 0; i < 3; i++){
            currentFace->halfedge[i].vertex->neighborHe = &(currentFace->halfedge[i]);
        }

        //If half edge collapsed has a mate i.e. not boundary then restore the face of mate as well
        if(hepCollapsed->mate != NULL){
            FaceIter oppositeFace = hepCollapsed->mate->face;

            oppositeFace->isActive = true;
            n_active_faces++;

            if(hepCollapsed->mate->next->mate != NULL) {
                hepCollapsed->mate->next->mate->mate = hepCollapsed->mate->next;
                hepCollapsed->mate->next->mate->edge = hepCollapsed->mate->next->edge;
            } 
            if(hepCollapsed->mate->prev->mate != NULL) {
                hepCollapsed->mate->prev->mate->mate = hepCollapsed->mate->prev;
                hepCollapsed->mate->prev->mate->edge = hepCollapsed->mate->prev->edge;
            } 

            hepCollapsed->mate->prev->edge->isActive = true;

            if(hepCollapsed->mate->prev->edge->halfedge[0] == hepCollapsed->mate->prev->mate) 
                hepCollapsed->mate->prev->edge->halfedge[1] = hepCollapsed->mate->prev;
            else {
                hepCollapsed->mate->prev->edge->halfedge[0] = hepCollapsed->mate->prev;
                /*if (hepCollapsed->mate->prev->edge->halfedge[1] != nullptr)
                {
                    hepCollapsed->mate->prev->edge->halfedge[1] = hepCollapsed->mate->prev->mate;
                }*/
            }

            if(hepCollapsed->mate->next->edge->halfedge[0] == hepCollapsed->mate->next->mate) 
                hepCollapsed->mate->next->edge->halfedge[1] = hepCollapsed->mate->next;
            else {
                hepCollapsed->mate->next->edge->halfedge[0] = hepCollapsed->mate->next;
                /*if (hepCollapsed->mate->prev->edge->halfedge[1] != nullptr)
                {
                    hepCollapsed->mate->prev->edge->halfedge[1] = hepCollapsed->mate->prev->mate;
                }*/
            }                  
                

            for(int i = 0; i < 3; i++){
                oppositeFace->halfedge[i].vertex->neighborHe = &(oppositeFace->halfedge[i]);
            }

        }

        //Restoring the halfedges for v0
        for(size_t i = 0; i < vertexSplitTarget.top().halfedgesAroundV0.size(); i++){
            //std::cout << "Does this run\n";
            HalfEdge *hep = vertexSplitTarget.top().halfedgesAroundV0[i];
            hep->vertex = v0;
        }


        //Restoring normal vectors
        UpdateNormalsAroundVertex(v0);
        UpdateNormalsAroundVertex(v1);

        //std::cerr << vertexSplitTarget.top().ei->halfedge[0]->edge->id << std::endl;
        //std::cerr << vertexSplitTarget.top().ei->halfedge[1]->edge->id << std::endl;
        //std::cerr << hepCollapsed->edge->id << std::endl;

        vertexSplitTarget.pop();
        
    } else {
        std::cerr << "Nothing to restore\n";
    }
}

//-----------------------------------------------------------------  

void Simplification::ControlLevelOfDetail(int step)
{
    int n_target_faces = mesh->n_faces*pow(0.95, step);
    /*int n_target_faces = mesh->n_faces;

    if (step >= 5)
    {
        n_target_faces -= 2;   
    }*/
    

    std::cerr << "step " << step << " " << n_target_faces << " " << mesh->n_faces << std::endl;

    if(n_target_faces < n_active_faces){
        while(n_target_faces < n_active_faces) {
            if(EdgeCollapse() == false) break;
            //std::cerr << n_active_faces << std::endl;
        } 
    }else if(n_target_faces > n_active_faces){
        while(n_target_faces > n_active_faces) {
            VertexSplit();
            //std::cerr << n_active_faces << std::endl;
        } 
    }

    std::cerr << "Current Active Faces: " << n_active_faces << std::endl;
    if (!validateMesh(this->mesh))
    {
        std::cerr << "Mesh is not valid!\n";
    }
    

}