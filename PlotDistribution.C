void PlotDistributions(TString filename){
  
  std::vector<TString> recon_method = {"Truth", "Electron", "JB", "DA", "Sigma", "ESigma"};

  // Open the file and retrieve the chain
  auto tree = new TChain("events");
  tree->Add(filename);
  
  for (auto method : recon_method){
    
    auto canvas = new TCanvas();
    canvas->Divide(2,2);
    
    TString branch_name;
    canvas->cd(1);
    
    // Draw a histogram for each variable as reconstructed by each method
    branch_name = TString::Format("InclusiveKinematics%s.x",method.Data());
    tree->Draw(branch_name);
    canvas->cd(2);
    branch_name = TString::Format("InclusiveKinematics%s.y",method.Data());
    tree->Draw(branch_name);
    canvas->cd(3);
    branch_name = TString::Format("InclusiveKinematics%s.Q2",method.Data());
    tree->Draw(branch_name);
    canvas->cd(4);
    branch_name = TString::Format("InclusiveKinematics%s.W",method.Data());
    tree->Draw(branch_name);

    branch_name = TString::Format("InclusiveKinematics%s.png",method.Data());
    canvas->Print(branch_name); // Write the canvases to a pdf file
  }
}