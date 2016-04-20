//Merge files using TFileMerger                                                                                                                                                                            
//Author: M. Verweij                                                                                                                                                                                       

void mergeFileMerger(std::vector<std::string> &infnames,  size_t startFile = 0, size_t nfiles = 0, const char *outname = "mergedForest.root") {

  if(nfiles>0) {
    nfiles = startFile + nfiles;
    printf("files2 %zu \n", nfiles);
    if(nfiles>infnames.size()) nfiles = infnames.size();
    printf("files2 %zu \n", nfiles);
  } else {
    startFile = 0;
    nfiles = infnames.size();
  }
  Printf("start file: %zu  nfiles: %zu",startFile,nfiles);

  //Create the TFileMerger instant                                                                                                                                                                         
  TFileMerger m(kFALSE);
  Int_t ma = 0;
  for(size_t i=startFile; i<nfiles; i++) {
    m.AddFile(infnames[i].c_str());
    ma++;
  }
  //m.SetNotrees(kTRUE);                                                                                                                                                                                   

  if(ma>0) {
    printf("outname %s \n", outname);
    m.OutputFile(outname);//"mergedForest.root");
    m.Merge();
  }
}
