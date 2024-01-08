/*
 * This file is the entry point into the reference implementation of Mint3D
 * Author: Josh Vekhter 
 */



#include "polyscope/polyscope.h"


// #include "VizHook.h"



////////////////////////////////////////////
////// A high level plan.  We Include a number of scripts for figures.  
////// Just turn one on at a time, and initize the Mint3DHook with the corresponding problem formulation.
////// There's a bit of boiler plate, sorry.  
////// Basically, according to the oracle, we store everything in an app state. 
////// The code is connected to an autodiff framework called TinyAD which supports with sparse hessians 
////// Each of these hooks sets up a specific problem, and ultimately corresponds to a different paper figure.
////// 
////// gpt said the app state looks like a mess, but I don't agree.: 
//////
//////  It's convenient. - ME  
//////  I think it's a good idea to have a single place to store all the data. -GPT 
//////  So you changed your mind? - ME
//////  No, I always thought it was a good idea. - GPT
//////
////// It could be organized better though, no? - ME 
////// I think it's organized fine. - GPT
//////
////////////////////////////////////////////


#include "MiNT_mesh.h"
// #include "MiNT_krushkal_rank2.h"

// #include "MiNT_krushkal_rank3.h"

// #include "MiNT_rank2_cylinder_example.h"



////////////////////////////////////////////


#include "Mint3DHook.h"

#include <thread>

#include "args.hxx"



static Mint3DHook *hook = NULL;


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

    std::cout << "try to reset" << std::endl;
    hook->reset();
    hook->appState->solveStatus = "reset state";



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
                // hook->appState->turn logging on 
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

  args::Flag headless(parser, "headless_mode", "This will avoid all of the polyscope calls and will immediately run", {'h', "headless"});
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

  if(!headless)
  {


  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  // polyscope::options::autocenterStructures = true;
// polyscope::options::autoscaleStructures = true;

  // Initialize polyscope

  // polyscope::init("openGL_mock");
    polyscope::init();

  }
  else
  {
    std::cout << "running MiNT 3D in headless mode" << std::endl;
  }



////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////

  // hook = static_cast<Mint3DHook*>(new MiNT_mesh<9>());
    hook = static_cast<Mint3DHook*>(new MiNT_mesh());

    // hook = static_cast<Mint3DHook*>(new MiNT_krushkal_rank2());
  // hook = static_cast<Mint3DHook*>(new MiNT_krushkal_rank3());
//
    // hook = static_cast<Mint3DHook*>(new MiNT_rank2_cylinder_example());

////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////

  hook->appState->headless_mode = headless;

  std::cout << "directoryPath: " << hook->appState->directoryPath << std::endl;
  hook->reset();
  std::cout << "nvars in opt: after reset" << hook->opt->get_num_vars() << std::endl; 

  hook->appState->directoryPath = args::get(inDir);

  if ( !( hook->appState->directoryPath.empty() )) {
    hook->appState->shouldReload = true;
    hook->appState->loadedPreviousRun = true;
  }
  else 
  {
    hook->reset();
  }

  if(!headless)
  {

  polyscope::state::userCallback = drawGUICallback;
  polyscope::options::programName = "gpgpt - MINT3D";
  polyscope::options::verbosity = 1;
  polyscope::options::transparencyRenderPasses = 4;
  polyscope::view::resetCameraToHomeView();
  polyscope::show();


  }
  else
  {
    hook->appState->headless_mode = true;
    // hook->run();

    while(hook->appState->keepSolving)
    {
      hook->simulateOneStep();
      hook->updateRenderGeometry();
    }

  }

  return 0;
}