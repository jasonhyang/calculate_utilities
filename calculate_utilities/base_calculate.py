from .calculate_smoothingFunctions import calculate_smoothingFunctions
from .calculate_statisticsSampledPoints import calculate_statisticsSampledPoints
from .calculate_statistics import calculate_statistics
from .calculate_biomass import calculate_biomass
from .calculate_clustering import calculate_clustering
from .calculate_curveFitting import calculate_curveFitting
from .calculate_histogram import calculate_histogram

class base_calculate(calculate_statistics,
                     calculate_smoothingFunctions,
                     calculate_statisticsSampledPoints,
                     calculate_biomass,
                     calculate_clustering,
                     calculate_curveFitting,
                     calculate_histogram):
    def __init__(self):
        self.data=[];
    # other
    def null(self, A, eps=1e-6):
        u, s, vh = numpy.linalg.svd(A,full_matrices=1,compute_uv=1)
        null_rows = [];
        rank = numpy.linalg.matrix_rank(A)
        for i in range(A.shape[1]):
            if i<rank:
                null_rows.append(False);
            else:
                null_rows.append(True);
        null_space = scipy.compress(null_rows, vh, axis=0)
        return null_space.T