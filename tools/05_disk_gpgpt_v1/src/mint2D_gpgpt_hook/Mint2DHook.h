#ifndef MINT2DHOOK_H
#define MINT2DHOOK_H

#include "AppState.h" // Include AppState definition
#include <Eigen/Dense>
#include <string>

class Mint2DHook {
public:
    Mint2DHook(AppState* state);
    virtual ~Mint2DHook();

    void drawGUI();
    void updateRenderGeometry();
    void renderRenderGeometry();
    void initSimulation();
    bool simulateOneStep();

private:
    AppState* appState; // Pointer to AppState instance

    void updateOptimizationParameters();
    void checkAndUpdateConvergence(double decrement, double energy);
    void updateVisualizationData();
    void finalizeIteration();

    // Other private member variables and functions as needed
};

#endif // MINT2DHOOK_H
