#include "GraphicsApp.hpp"

int main() {
    
    GraphicsApp app(640, 480, -0.5f, 0.0f);

    app.Initialize();
    app.MainLoop();
    app.CleanUp();

    return 0;
}
