import csv
import numpy

def generate_row(numCols):
  rtn = [0.0] * numCols
  for i in range(numCols):
    rtn[i] = numpy.random.uniform(-999999, 999999)
  return rtn

class Main:
  filename = "trainingData.csv"
  f = open(filename, 'wb')
  writer = csv.writer(f)
  numRows = 32000000
  numCols = 8
  for i in range(numRows):
    row = generate_row(numCols);
    writer.writerow(row)
