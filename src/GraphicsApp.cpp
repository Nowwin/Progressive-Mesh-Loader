#include "GraphicsApp.hpp"


int GNumberOfVertices = 0;
bool wireframeMode = false;

Mesh bunny = Mesh();

// Camera
Camera gCamera;

void GraphicsApp::GetOpenGLVersionInfo() {
    std::cout << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "Shading Language: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
}

void GraphicsApp::IntializeProgram() {

    std::cout << std::filesystem::current_path().string() << std::endl; 

    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        std::cout << "SDL did not initiliaze" << std::endl;
        exit(1);
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    gGraphicsApplicationWindow = SDL_CreateWindow(
        "SDL2 Window",               // window title
        0,      // x position, centered
        0,      // y position, centered
        m_screenWidth,                         // width, in pixels
        m_screenHeight,                         // height, in pixels
        SDL_WINDOW_OPENGL             // flags
    );

    if (gGraphicsApplicationWindow == nullptr)
    {
        std::cout << "Window did not setup" << std::endl;
        exit(1);
    }

    gOpenGlContext = SDL_GL_CreateContext(gGraphicsApplicationWindow);

    if (gOpenGlContext == nullptr)
    {
        std::cout << "Graphics did not load" << std::endl;
        exit(1);
    }
    
    if (!gladLoadGLLoader(SDL_GL_GetProcAddress))
    {
        std::cout << "glad was not initialized" << std::endl;
        exit(1);
    }

    GetOpenGLVersionInfo();
    
}

void GraphicsApp::VertexSpecification() {
    
    char filename[] = "./src/ObjFiles/ver1.obj";
    bunny.ConstructMeshDataStructure(filename);
    Simplification simple = Simplification();
    simple.InitSimplification(&bunny);
    simple.ControlLevelOfDetail(50);
    simple.ControlLevelOfDetail(40);
    Mesh* bunnyModified = simple.GetModifiedMesh();
    const std::vector<GLfloat> vertexData(bunnyModified->GetVertexData());

    glGenVertexArrays(1, &gVertexArrayObject);
    //Selecting
    glBindVertexArray(gVertexArrayObject);

    glGenBuffers(1, &gVertexBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexBufferObject);
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(GLfloat), vertexData.data(), GL_STATIC_DRAW);

   

    const std::vector<GLuint> indexBufferData(bunnyModified->GetIndexData());

    GNumberOfVertices = indexBufferData.size();
                                                
    glGenBuffers(1, &gIndexBufferObject);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBufferObject);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexBufferData.size() * sizeof(GLuint), indexBufferData.data(), GL_STATIC_DRAW);


    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT) * 6, (void*)0);
    
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT) * 6, (GLvoid*)(sizeof(GL_FLOAT)*3));

    glBindVertexArray(0);    
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    
}

void GraphicsApp::LoadShaders() {
    std::string vertexShaderSource = ShaderLoader::LoadShaderAsString(std::filesystem::current_path().string() + "/shaders/vert.glsl");
    std::string fragmentShaderSource = ShaderLoader::LoadShaderAsString(std::filesystem::current_path().string() + "/shaders/frag.glsl");
    gGraphicsPipelineShaderProgram = ShaderLoader::CreateShaderProgram(vertexShaderSource, fragmentShaderSource);
}


void GraphicsApp::Input(){
    // Two static variables to hold the mouse position
    static int mouseX=m_screenWidth/2;
    static int mouseY=m_screenHeight/2; 

    SDL_Event event;

    // (1) Handle Input
    // Start our event loop
    while(SDL_PollEvent(&event)){
        // Handle each specific event
        if(event.type == SDL_QUIT){
            gQuit= true;
        }

        if(event.type == SDL_KEYDOWN && (event.key.keysym.sym == SDLK_ESCAPE || event.key.keysym.sym == SDLK_q)){
			std::cout << "ESC: Goodbye! (Leaving MainApplicationLoop())" << std::endl;
            gQuit = true;
        }

        if(event.type==SDL_MOUSEMOTION){
            // Capture the change in the mouse position
            mouseX+=event.motion.xrel;
            mouseY+=event.motion.yrel;
            gCamera.MouseLook(mouseX,mouseY);
        }

        const Uint8 *state = SDL_GetKeyboardState(NULL);
        if (state[SDL_SCANCODE_UP]) {
            gCamera.MoveForward(0.1f);
        }
        if (state[SDL_SCANCODE_DOWN]) {
            gCamera.MoveBackward(0.1f);
        }
        if (state[SDL_SCANCODE_LEFT]) {
            gCamera.MoveLeft(0.1f);
        }
        if (state[SDL_SCANCODE_RIGHT]) {
            gCamera.MoveRight(0.1f);
        }

        if (event.type == SDL_KEYDOWN)
        {
            switch (event.key.keysym.sym) {
                    case SDLK_w:
                        // Handle 'W' key press here
                        printf("W key pressed!\n");
                        wireframeMode = !wireframeMode;
                        break;
                    default:
                        break;
            
             }
        }
        
    }
}

void GraphicsApp::PreDraw() {
    //Setting opengl state
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glViewport(0, 0, m_screenWidth, m_screenHeight);
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(gGraphicsPipelineShaderProgram);

    
    // Model transformation by translating our object into world space
    glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));

    GLint u_modelMatrixLocation = glGetUniformLocation(gGraphicsPipelineShaderProgram, "u_ModelMatrix");
    if(u_modelMatrixLocation >=0){
        glUniformMatrix4fv(u_modelMatrixLocation,1,GL_FALSE,&model[0][0]);
    }else{
        std::cout << "Could not find u_ModelMatrix, maybe a mispelling?\n";
        exit(EXIT_FAILURE);
    }

    // Update the View Matrix
    GLint u_ViewMatrixLocation = glGetUniformLocation(gGraphicsPipelineShaderProgram,"u_ViewMatrix");
    if(u_ViewMatrixLocation>=0){
        glm::mat4 viewMatrix = gCamera.GetViewMatrix();
        glUniformMatrix4fv(u_ViewMatrixLocation,1,GL_FALSE,&viewMatrix[0][0]);
    }else{
        std::cout << "Could not find u_ModelMatrix, maybe a mispelling?\n";
        exit(EXIT_FAILURE);
    }

    //Perspective

    glm::mat4 perspective = glm::perspective(glm::radians(45.0f), (float)m_screenWidth/(float)m_screenHeight, 0.1f, 20.0f);

    GLint u_perspectiveLocation = glGetUniformLocation(gGraphicsPipelineShaderProgram, "u_PerspectiveMatrix");

    if (u_perspectiveLocation >= 0)
    {
        glUniformMatrix4fv(u_perspectiveLocation, 1, GL_FALSE, &perspective[0][0]);
    } else {
        std::cout << "No location found" << std::endl;
    }

    GLint isWireframeModeLocation = glGetUniformLocation(gGraphicsPipelineShaderProgram, "isWireframeMode");

    //WireFrame Mode

    if (isWireframeModeLocation >= 0)
    {
        glUniform1i(isWireframeModeLocation, wireframeMode);
    } else {
        std::cout << "No location found" << std::endl;
    }
    
}

void GraphicsApp::Draw() {

    glBindVertexArray(gVertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexBufferObject);

    //glDrawArrays(GL_TRIANGLES, 0, 6);
    if (wireframeMode)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    
    

    glDrawElements(GL_TRIANGLES, GNumberOfVertices, GL_UNSIGNED_INT, 0);
    

}

void GraphicsApp::Initialize() {
    IntializeProgram();
    
    VertexSpecification();
    
    LoadShaders();
}

void GraphicsApp::MainLoop() {

    SDL_WarpMouseInWindow(gGraphicsApplicationWindow,m_screenWidth/2,m_screenHeight/2);
    SDL_SetRelativeMouseMode(SDL_TRUE);

    while (!gQuit) {
        Input();
        
        PreDraw();

        Draw();

        SDL_GL_SwapWindow(gGraphicsApplicationWindow);
    }
}

void GraphicsApp::CleanUp() {
    SDL_DestroyWindow(gGraphicsApplicationWindow);
    SDL_Quit();
}

