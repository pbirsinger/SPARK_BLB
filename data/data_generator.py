import csv
import numpy
import pickle


def generate_row(numCols):
  rtn = [0.0] * numCols
  for i in range(numCols):
    rtn[i] = numpy.random.uniform(-999999, 999999)
  return rtn

class Main:
  """
  filename = "eight_dim_k_means_data.csv"
  f = open(filename, 'wb')
  writer = csv.writer(f)
  numRows = 100000
  numCols = 8
  for i in range(numRows):
    row = generate_row(numCols);
    writer.writerow(row)
    """
  M = numpy.random.rand(25 * 10 ** 4, 1)
  outfile = open('onegig.pkl', 'wb')
  pickle.dump(M, outfile)
  outfile.close();

  infile = open('onegig.pkl', 'rb')
  N = pickle.load(infile)
  infile.close()

  print len(N)
