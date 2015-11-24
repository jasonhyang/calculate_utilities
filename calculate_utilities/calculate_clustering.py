from .calculate_dependencies import *
class calculate_clustering():
    def __init__(self):
        self.data=[];
    # heatmap
    def calculate_heatmap(self,data_I,row_labels_I,column_labels_I,
                row_pdist_metric_I='euclidean',row_linkage_method_I='complete',
                col_pdist_metric_I='euclidean',col_linkage_method_I='complete'):
        '''Generate a heatmap using pandas and scipy'''

        """dendrogram documentation:
        linkage Methods:
        single(y)	Performs single/min/nearest linkage on the condensed distance matrix y
        complete(y)	Performs complete/max/farthest point linkage on a condensed distance matrix
        average(y)	Performs average/UPGMA linkage on a condensed distance matrix
        weighted(y)	Performs weighted/WPGMA linkage on the condensed distance matrix.
        centroid(y)	Performs centroid/UPGMC linkage.
        median(y)	Performs median/WPGMC linkage.
        ward(y)	Performs Ward's linkage on a condensed or redundant distance matrix.
        Output:
        'color_list': A list of color names. The k?th element represents the color of the k?th link.
        'icoord' and 'dcoord':  Each of them is a list of lists. Let icoord = [I1, I2, ..., Ip] where Ik = [xk1, xk2, xk3, xk4] and dcoord = [D1, D2, ..., Dp] where Dk = [yk1, yk2, yk3, yk4], then the k?th link painted is (xk1, yk1) - (xk2, yk2) - (xk3, yk3) - (xk4, yk4).
        'ivl':  A list of labels corresponding to the leaf nodes.
        'leaves': For each i, H[i] == j, cluster node j appears in position i in the left-to-right traversal of the leaves, where \(j < 2n-1\) and \(i < n\). If j is less than n, the i-th leaf node corresponds to an original observation. Otherwise, it corresponds to a non-singleton cluster."""

        #parse input into col_labels and row_labels
        #TODO: pandas is not needed for this.
        mets_data = pd.DataFrame(data=data_I, index=row_labels_I, columns=column_labels_I)

        mets_data = mets_data.dropna(how='all').fillna(0.)
        #mets_data = mets_data.replace([np.inf], 10.)
        #mets_data = mets_data.replace([-np.inf], -10.)
        col_labels = list(mets_data.columns)
        row_labels = list(mets_data.index)

        #perform the custering on the both the rows and columns
        dm = mets_data
        D1 = squareform(pdist(dm, metric=row_pdist_metric_I))
        D2 = squareform(pdist(dm.T, metric=col_pdist_metric_I))

        Y = linkage(D1, method=row_linkage_method_I)
        Z1 = dendrogram(Y, labels=dm.index)

        Y = linkage(D2, method=col_linkage_method_I)
        Z2 = dendrogram(Y, labels=dm.columns)

        #parse the output
        col_leaves = Z2['leaves'] # no hclustering; same as heatmap_data['col']
        row_leaves = Z1['leaves'] # no hclustering; same as heatmap_data['row']
        col_colors = Z2['color_list'] # no hclustering; same as heatmap_data['col']
        row_colors = Z1['color_list'] # no hclustering; same as heatmap_data['row']
        col_icoord = Z2['icoord'] # no hclustering; same as heatmap_data['col']
        row_icoord = Z1['icoord'] # no hclustering; same as heatmap_data['row']
        col_dcoord = Z2['dcoord'] # no hclustering; same as heatmap_data['col']
        row_dcoord = Z1['dcoord'] # no hclustering; same as heatmap_data['row']
        col_ivl = Z2['ivl'] # no hclustering; same as heatmap_data['col']
        row_ivl = Z1['ivl'] # no hclustering; same as heatmap_data['row']

        #heatmap data matrix
        heatmap_data_O = []
        for i,r in enumerate(mets_data.index):
            for j,c in enumerate(mets_data.columns):
                #heatmap_data.append({"col": j+1, "row": i+1, "value": mets_data.ix[r][c]})
                tmp = {"col_index": j, "row_index": i, "value": mets_data.ix[r][c],
                       'row_label': r,'col_label':c,
                       'col_leaves':col_leaves[j],
                        'row_leaves':row_leaves[i],
                        'col_pdist_metric':col_pdist_metric_I,
                        'row_pdist_metric':row_pdist_metric_I,
                        'col_linkage_method':col_linkage_method_I,
                        'row_linkage_method':row_linkage_method_I,
                       };
                heatmap_data_O.append(tmp)

        dendrogram_col_O = {'leaves':col_leaves,
                        'icoord':col_icoord,
                        'dcoord':col_dcoord,
                        'ivl':col_ivl,
                        'colors':col_colors,
                        'pdist_metric':col_pdist_metric_I,
                        'linkage_method':col_linkage_method_I};

        dendrogram_row_O = {
                        'leaves':row_leaves,
                        'icoord':row_icoord,
                        'dcoord':row_dcoord,
                        'ivl':row_ivl,
                        'colors':row_colors,
                        'pdist_metric':row_pdist_metric_I,
                        'linkage_method':row_linkage_method_I};
        return heatmap_data_O,dendrogram_col_O,dendrogram_row_O;
    
    def calculate_clusters(self,data_I,cluster_method_I):
        '''call a clustering algorithm on the data
        Cluster methods available:
        sklearn.cluster.KMeans: 
        The simplest, yet effective clustering algorithm. Needs to be provided with the number of clusters in advance, and assumes that the data is normalized as input (but use a PCA model as preprocessor).
        sklearn.cluster.MeanShift: 
        Can find better looking clusters than KMeans but is not scalable to high number of samples.
        sklearn.cluster.DBSCAN: 
        Can detect irregularly shaped clusters based on density, i.e. sparse regions in the input space are likely to become inter-cluster boundaries. Can also detect outliers (samples that are not part of a cluster).
        sklearn.cluster.AffinityPropagation: 
        Clustering algorithm based on message passing between data points.
        sklearn.cluster.SpectralClustering: 
        KMeans applied to a projection of the normalized graph Laplacian: finds normalized graph cuts if the affinity matrix is interpreted as an adjacency matrix of a graph.
        sklearn.cluster.Ward: 
        Ward implements hierarchical clustering based on the Ward algorithm, a variance-minimizing approach. At each step, it minimizes the sum of squared differences within all clusters (inertia criterion).
        '''
        pass;

    #k-means
    def calculate_kmeans(self):
        '''calculate the kmeans clusters
        sklearn.cluster.k_means(X, n_clusters, init='k-means++', precompute_distances='auto', n_init=10, max_iter=300, verbose=False, tol=0.0001, random_state=None, copy_x=True, n_jobs=1, return_n_iter=False)
        
        E.g.
        from sklearn.datasets import make_blobs
        X, y = make_blobs(random_state=42)
        X.shape

        plt.scatter(X[:, 0], X[:, 1])

        from sklearn.cluster import KMeans

        kmeans = KMeans(n_clusters=3, random_state=42)
        labels = kmeans.fit_predict(X)

        plt.scatter(X[:, 0], X[:, 1], c=labels)
        '''
        pass

    def score_clusters(self,data_I,labels_I,metric_I):
        '''Score the clustering algorithm
        Methods:
        from sklearn.metrics import confusion_matrix,
            accuracy_score
            adjusted_rand_score
        
        '''

    #SOM