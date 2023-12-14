/*
 * Adapted from:
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */



#include "polyscope/polyscope.h"

// #include "VizHook.h"

#include "Solve_L2_newton_rank1.h"
#include "Mint2DHook.h"

#include <thread>

#include "args.hxx"



static Mint2DHook *hook = NULL;


void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
    {
        hook->run();
    }
    else
    {
        hook->pause();
    }
        
}

void resetSimulation()
{
    
    if (!hook)
        return;

  // double w_bound_prev = hook->w_bound;
  // double w_smooth_prev = hook->w_smooth;
  // double w_curl_prev = hook->w_curl;
  // double w_s_perp_prev = hook->w_s_perp;

    std::cout << "try to reset" << std::endl;
    hook->reset();
    hook->appState->solveStatus = "reset state";

  // hook->w_bound = w_bound_prev;
  // hook->w_smooth = w_smooth_prev;
  // hook->w_curl = w_curl_prev;
  // hook->w_s_perp = w_s_perp_prev;

}

void drawGUICallback()
{
	ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
							   // instead of full width. Must have 
							   // matching PopItemWidth() below.

    if(hook->showSimButtons() || true)
    {
        if (ImGui::CollapsingHeader("Start Simulation.", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Run/Pause Sim"))
            {
                toggleSimulation();
            }
            if (ImGui::Button("Reset Sim"))
            {
                resetSimulation();
            }
            if (ImGui::Button("Take one step"))
            {
                if (hook->isPaused())
                {
                  hook->simulateOneStep();
                  hook->updateRenderGeometry();
                }
                  
            }
        }
    }
    hook->drawGUI();
	ImGui::PopItemWidth();
    hook->render();
}

int main(int argc, char **argv) {
  // Configure the argument parser

// Comment out parser.  Should be easy to add in commandline scripting with this though.

  args::ArgumentParser parser(("This is an example optimization that studies integrability in 2d.\n "
                              "Built on top of an LLM assisted codebase called gpgpt. \n\n\n"
                              ""));
  args::Positional<std::string> inDir(parser, "out_dir", "load a directory to restart optimization");
  args::Positional<std::string> inObj(parser, "in_obj", "load a surface mesh to solve if out_dir is empty");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
    // std::cout << parser;
  } catch (args::Help&) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError& e) {
    std::cerr << e.what() << std::endl;

    std::cerr << parser;
    return 1;
  }



  


  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  polyscope::options::autocenterStructures = true;
// polyscope::options::autoscaleStructures = true;

  // Initialize polyscope
  polyscope::init();

  hook = static_cast<Mint2DHook*>(new Solve_L2_newton_rank1());



  std::cout << "directoryPath: " << hook->appState->directoryPath << std::endl;
  hook->reset();
  std::cout << "nvars in opt: after reset" << hook->opt->get_num_vars() << std::endl; 

  hook->appState->directoryPath = args::get(inDir);

  if ( !( hook->appState->directoryPath.empty() )) {
    // hook->initSimulation();
    // hook->appState->objFilePath = args::get(inObj);
    hook->appState->shouldReload = true;
  }

  // hook->initSimulation();
  // std::cout << "nvars in opt: outside after init" << hook->opt->get_num_vars() << std::endl; 

  polyscope::state::userCallback = drawGUICallback;
  polyscope::options::programName = "gpgpt - MINT2D";
  polyscope::options::verbosity = 1;

//   polyscope::GroundPlaneMode::None;
// polyscope::options::groundPlaneEnabled = false;

// polyscope::screenshotExtension = ".jpg";


  polyscope::show();

  return 0;
}