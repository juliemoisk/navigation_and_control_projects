import numpy as np

class KMeans(object):
    def init(self, K, init):
        self.K = K                    
        self.centroids = np.array(init) 
    
    def fit(self, X):
        while True:
            distances = np.linalg.norm(X[:, np.newaxis] - self.centroids, axis=2)
            
            labels = np.argmin(distances, axis=1)
            
            new_centroids = np.array([
                X[labels == i].mean(axis=0) if np.any(labels == i) else self.centroids[i]
                for i in range(self.K)
            ])
            
            shifts = np.linalg.norm(new_centroids - self.centroids, axis=1)
            if np.all(shifts <= 0.001):
                break
            
            self.centroids = new_centroids
    
    def predict(self, X):
        distances = np.linalg.norm(X[:, np.newaxis] - self.centroids, axis=2)
        labels = np.argmin(distances, axis=1)
        return np.array(labels)