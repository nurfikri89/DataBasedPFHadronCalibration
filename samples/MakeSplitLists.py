import itertools
import os
import glob

# def chunks(l, n):
#   """Yield successive n-sized chunks from l."""
#   for i in range(0, len(l), n):
#     yield l[i:i + n]

# def main(txtFileName, nGroupFiles=10):
#   print(f"Split for {txtFileName}")
#   listOfTxtFiles = glob.glob(txtFileName+".txt")
#   for txtFile in listOfTxtFiles:
#     sampleName = txtFile.replace(".txt","")
#     with open(txtFile) as f:
#       files = f.read().splitlines()
#       filesChunks = list(chunks(files,nGroupFiles))
#       i=0
#       for fileChunk in filesChunks:
#         txtFileSplitName = sampleName+"_part"+str(i)+".txt"
#         fOut = open(txtFileSplitName, "w")
#         for line in fileChunk:
#           # write line to output file
#           fOut.write(line + "\n")
#         fOut.close()
#         i += 1

def main(sampleName, nPartFiles=2):
  txtFileName=f"{sampleName}.txt"
  print(f"Split for {txtFileName}")
  with open(txtFileName, "r", encoding="utf-8") as f:
    lines = [line.strip() for line in f if line.strip()]  # remove blanks

  total_lines = len(lines)
  if total_lines == 0:
    print("The input file is empty.")
    return

  # Calculate how many lines per chunk
  chunk_size = total_lines // nPartFiles
  remainder = total_lines % nPartFiles

  start = 0
  for i in range(nPartFiles):
    # Distribute remainder (extra lines) across the first few chunks
    end = start + chunk_size + (1 if i < remainder else 0)
    chunk = lines[start:end]
    start = end

    # Write chunk to output file
    output_filename = f"{sampleName}_part{i}.txt"
    with open(output_filename, "w", encoding="utf-8") as out:
      out.write("\n".join(chunk))
    print(f"Created: {output_filename} ({len(chunk)} lines)")

if __name__ == "__main__":

  txtFiles=[
    ("MC24_SingleNeutrino",5),
  ]

  for txtFile in txtFiles:
    main(txtFile[0],txtFile[1])

