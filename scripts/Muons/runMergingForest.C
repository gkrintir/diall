#include "mergeFileMerger.C"
#include <sstream>
#include <string>
#include <vector>
#include "TString.h"
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atol */

std::vector<std::string> CollectFiles(TString list) {

  std::vector<std::string> urls;

  ifstream inputFile (list.Data());
  string line;
  if (!inputFile.is_open())
    cout<<"Could not open file\n";
  else 
  {
    do
    {
      getline(inputFile,line); 
      if (inputFile)
      {
	stringstream s(line); 
	if (s.good())
	{
	  printf("line: %s \n", line.data());
	  urls.push_back(line.data());
	}
      }
    }while(inputFile);
  }
  inputFile.close();
  return urls;

}

void runMergingForest(TString list = "output.list", size_t startFile = 0, size_t nfiles = 0, const char *outname = "mergedForest.root") {

  std::vector<std::string> urls = CollectFiles(list);

  // for (std::vector<std::string>::const_iterator i = urls.begin(); i != urls.end(); ++i)                                                                                                                 
  //   std::cout << *i << std::endl;                                                                                                                                                                       


  // TString macro = Form(merge.C);                                                                                                                                                                       
  printf("files!! %zu \n", nfiles);
  mergeFileMerger(urls,startFile,nfiles,outname);

}
