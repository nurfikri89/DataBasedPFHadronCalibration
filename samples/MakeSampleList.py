import os
import subprocess
import sys, json

def main(datasetName,sampleName, USERdataset=False):

  list_file_path = []

  if type(datasetName) == str:
    if not(USERdataset):
      out = subprocess.check_output("dasgoclient -json --query 'file dataset=%s'" %(datasetName), shell=True)
    else:
      out = subprocess.check_output("dasgoclient -json --query 'file dataset=%s instance=prod/phys03'" %(datasetName), shell=True)

    list_jsonDict = json.loads(out) # Get list of dictionary. Each element of list is for a file.
    for data in list_jsonDict:
      file_path = data["file"][0]["name"]
      list_file_path.append(file_path)
  elif type(datasetName) == list:
    print(f"Have extension datasets for sample: ({sampleName})")
    for thisDataset in datasetName:
      print(thisDataset)
      if not(USERdataset):
        out = subprocess.check_output("dasgoclient -json --query 'file dataset=%s'" %(thisDataset), shell=True)
      else:
        out = subprocess.check_output("dasgoclient -json --query 'file dataset=%s instance=prod/phys03'" %(thisDataset), shell=True)

      list_jsonDict = json.loads(out) # Get list of dictionary. Each element of list is for a file.

      for data in list_jsonDict:
        file_path = data["file"][0]["name"]
        list_file_path.append(file_path)
  #
  #
  #
  fout = f"{sampleName}.txt"
  fo = open(fout, "w")
  print(f"Making {fout} with nfiles={len(list_file_path)}")
  for file_path in list_file_path:
    fo.write(file_path+'\n')
  fo.close()

datasetDict = {
"MC24_SingleNeutrino":"/SingleNeutrino_Par-E-10_gun/nbinnorj-RunIII2024Summer24NanoAODv15_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24C_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024C_NanoV15_v1_2024CDEReprocessing_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24D_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024D_NanoV15_v1_2024CDEReprocessing_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24E_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024E_NanoV15_v1_2024CDEReprocessing_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24F_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024F_NanoV15_v1_Prompt_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24G_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024G_NanoV15_v1_Prompt_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24H_ZeroBias"  :"/ZeroBias/nbinnorj-Run2024H_NanoV15_v1_Prompt_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24Iv1_ZeroBias":"/ZeroBias/nbinnorj-Run2024I_NanoV15_v1_Prompt_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
"Data24Iv2_ZeroBias":"/ZeroBias/nbinnorj-Run2024I_NanoV15_v2_Prompt_AODToPFNANO_v15p0-00000000000000000000000000000000/USER",
}


for key in datasetDict:
  datasetName=datasetDict[key]
  sampleName=key
  main(datasetName,sampleName,USERdataset=True)
  # main(datasetName,sampleName,USERdataset=False)


