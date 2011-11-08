import csv
import math
from random import randint

def recalc_mean(originalMean, originalNumItems, item):
  result = originalMean
  for i in range(len(result)):
    result[i] = result[i] * originalNumItems
    result[i] += item[i]
    result[i] = result[i] / (originalNumItems + 1)
  return result

def dist_to_mean(observation, mean):
  if (len(observation) != len(mean)):
    return 999999
  dist = 0.0
  for i in range(len(observation)):
    dist += (observation[i] - mean[i]) * (observation[i] - mean[i])
  dist = math.sqrt(dist)
  return dist

def k_means(numClusters, dimensions):

  filename = "eight_dim_k_means_data.csv"
  f = open(filename)
  observations = csv.reader(f)

  clusters = []
  #Create k clusters
  for i in xrange(numClusters):
    clusters.append(Cluster(dimensions))
  #Randomly assign each observation to a Cluster
  for observation in observations:
    for i in xrange(len(observation)):
      observation[i] = float(observation[i])
    clusters[randint(0, numClusters - 1)].add(observation)

  changed = True
  while (changed):
    changed = False
    for cluster1 in clusters:
      for observation in cluster1.get_items():
        distToCurrentClusterMean = 999999
        currentCluster = None
        bestCluster = None
        distToBestClusterMean = float("inf")
        for cluster in clusters:
          #Get distance to each cluster's mean
          dist = dist_to_mean(observation, cluster.get_mean())
          if (cluster.contains(observation)):
            distToCurrentClusterMean = dist
            currentCluster = cluster
          if dist <= distToBestClusterMean:
            distToBestClusterMean = dist
            bestCluster = cluster
        if distToBestClusterMean < distToCurrentClusterMean:
          #Move to bestCluster
          changed = True
          currentCluster.remove(observation)
          bestCluster.add(observation)
  return clusters

class Cluster():
  def __init__(self, dimensions):
    self.mean = [0.0] * dimensions
    self.items = []
    self.numItems = 0

  def add(self, item):
    self.numItems += 1
    self.items.append(item)
    self.mean = recalc_mean(self.mean, self.numItems - 1, item)

  def remove(self, item):
    self.numItems -= 1
    self.items.remove(item)
    self.mean = recalc_mean(self.mean, self.numItems - 1, item)

  def get_items(self):
    return self.items

  def get_mean(self):
    return self.mean

  def contains(self, item):
    return (item in self.items)

  def get_num_items(self):
    return self.numItems

class Main():
  numClusters = 3
  dimensions = 8
  clusters = k_means(numClusters, dimensions)
  for cluster in clusters:
    print cluster.get_num_items()
