import numpy as np

class sample(object):
  def __init__(self, X, n_subspace):
    self.idx_subspace = self.random_subspace(X, n_subspace)

  def __call__(self, X, y):
    idx_obj = self.bootstrap_sample(X)
    X_sampled, y_sampled = self.get_subsample(X, y, self.idx_subspace, idx_obj)
    return X_sampled, y_sampled

  @staticmethod
  def bootstrap_sample(X, random_state=42):
    return np.unique(np.random.choice(X.shape[0], X.shape[0]))

  @staticmethod
  def random_subspace(X, n_subspace, random_state=42):
    return np.sort(np.random.choice(X.shape[1], n_subspace, replace=False))

  @staticmethod
  def get_subsample(X, y, idx_subspace, idx_obj):
    return X[idx_obj][:, idx_subspace], y[idx_obj]
######################################################################################## первое
  import numpy as np
  from sklearn.datasets import load_iris
  from sklearn.tree import DecisionTreeClassifier
  from sklearn.metrics import accuracy_score
  from sklearn.model_selection import train_test_split
  from collections import Counter

  class sample:
    def __init__(self, X: np.array, n_subspace: int):
      self.idx_subspace = self.random_subspace(X, n_subspace)

    def __call__(self, X: np.array, y: np.array):
      idx_obj = self.bootstrap_sample(X)
      X_sampled, y_sampled = self.get_subsample(
        X, y, self.idx_subspace, idx_obj
      )
      return X_sampled, y_sampled

    @staticmethod
    def bootstrap_sample(X: np.array):
      return np.unique(np.random.choice(X.shape[0], X.shape[0]))

    @staticmethod
    def random_subspace(X: np.array, n_subspace: int):
      return np.sort(np.random.choice(X.shape[1], n_subspace, replace=False))

    @staticmethod
    def get_subsample(
            X: np.array, y: np.array, idx_subspace: np.array, idx_obj: np.array
    ):
      return X[idx_obj][:, idx_subspace], y[idx_obj]

  class random_forest:
    def __init__(
            self,
            n_estimators: int,
            max_depth: int,
            subspaces_dim: int,
            random_state: int,
    ):
      self.n_estimators = n_estimators
      self.max_depth = max_depth
      self.subspaces_dim = subspaces_dim
      self.random_state = random_state
      self.Classifier = []
      self.subspace_idx = []

    def fit(self, X: np.array, y: np.array):
      for _ in range(self.n_estimators):
        s = sample(X, self.subspaces_dim)
        X_sample, y_sample = s(X, y)
        self.Classifier.append(
          DecisionTreeClassifier(max_depth=self.max_depth).fit(
            X_sample, y_sample
          )
        )
        self.subspace_idx.append(s.idx_subspace)

    def predict(self, X: np.array):
      tmp = np.empty((self.n_estimators, X.shape[0]), dtype=np.float64)

      for i in range(self.n_estimators):
        tmp[i, :] = self.Classifier[i].predict(X[:, self.subspace_idx[i]])

      out = tmp.T
      predicted = np.zeros(X.shape[0], dtype=int)

      for i, x in enumerate(out):
        counter = Counter(x)
        predicted[i] = counter.most_common(1)[0][0]

      return predicted