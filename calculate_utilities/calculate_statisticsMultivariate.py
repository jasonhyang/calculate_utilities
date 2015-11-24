from .calculate_dependencies import *
class calculate_statisticsMultivariate():
    def __init__(self):
        self.data=[];

    def calculate_pca(self):
        '''calculate PCA
        sklearn.decomposition.PCA(n_components=None, copy=True, whiten=False)

        from sklearn.datasets import load_iris

        iris = load_iris()
        X, y = iris.data, iris.target

        from sklearn.decomposition import PCA
        pca = PCA()

        pca.fit(X_blob)

        X_pca = pca.transform(X_blob)

        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y, linewidths=0, s=30)
        plt.xlabel("first principal component")
        plt.ylabel("second principal component")

        '''
        pass;

    def calculate_manifold(self,data_I,method_I):
        '''Calculate the dimensionality reduction using manifold methods
        Methods:
        from sklearn.manifold import isomap, TSNE, ...
        '''

    def calculate_plsda(self):
        '''calculate PLS-DA using PLS
        '''
        pass;

    def calculate_pls(self):
        '''calculate PLS'''
        #break into training and testing data splits
        from sklearn.cross_validation import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

        #score the model

        #validate the model
        #cross validation the model
        #permute the model
        #make the final model
        pass;

    def calculate_pls_explainedVariance(self):
        '''calculate the explained variance of the pls scores'''
        pass;

    def calculate_pls_vip(self):
        '''calculate the variance of the pls scores and loadings'''
        pass;

    def crossValidate_pls(self,X,y,classifier,cv):
        '''cross validate the model
        METHOD:
        cross_val_score uses StratifiedKFold by default
        other cross_validation classes include KFold, StratifiedKFold, ShuffleSplit, LeavePLabelOut
        INPUT:
        E.g.
        from sklearn.datasets import load_iris
        from sklearn.neighbors import KNeighborsClassifier

        iris = load_iris()
        X, y = iris.data, iris.target

        classifier = KNeighborsClassifier()

        X = samples, array of shape len(n_samples),len(n_features)
        y = target, array of shape len(n_samples)
        classifier = class, classifier used
        cv = number of folds or object of class cross_validation

        where cv
        cv = KFold(n=len(x), shuffle=True)
        cv = StratifiedKFold(iris.target, n_folds=5)
        cv = ShuffleSplit(len(iris.target), n_iter=5, test_size=.2)
        ...
        
        '''
        

        from sklearn.cross_validation import cross_val_score
        scores = cross_val_score(classifier, X, y,cv=cv);

        pass;

    def gridSearch_parameters(self):
        '''perform a grid serch of the parameter space
        E.g.:
        from sklearn.grid_search import GridSearchCV
        from sklearn.svm import SVR
        param_grid = {'C': [0.001, 0.01, 0.1, 1, 10], 'gamma': [0.001, 0.01, 0.1, 1]}

        grid = GridSearchCV(SVR(), param_grid=param_grid, cv=cv, verbose=3)
        grid.fit(X, y)
        grid.predict(X)
        print(grid.best_score_)
        print(grid.best_params_)

        if over fitting is an issue
        from sklearn.cross_validation import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

        param_grid = {'C': [0.001, 0.01, 0.1, 1, 10], 'gamma': [0.001, 0.01, 0.1, 1]}
        cv = KFold(n=len(X_train), n_folds=10, shuffle=True)

        grid = GridSearchCV(SVR(), param_grid=param_grid, cv=cv)

        grid.fit(X_train, y_train)
        grid.score(X_test, y_test)
        '''

    def permute_pls(self):
        '''perform a permutation test on a cross-validated model to test for significance

        from sklearn.cross_validation import permutation_test_score
        
        '''
        pass;

    def calculate_svd(self):
        '''calculate SVD'''
        pass;

    def convert_listDict2dataMatrix(self,data_I,row_label_I,column_label_I,
                                    row_variables_I=[],column_variables_I=[]
                                    ):
        '''convert a list of dictionary rows to a numpy array
        INPUT:
        data_I = [{}]
        row_label_I = column_id of the row labels
        column_label_I = column_id of the column labels

        OPTIONAL INPUT:
        row_variables_I = list of keys to extract out with the rows
        column_variables_I = list of keys to extract out with the columns

        OUTPUT:
        data_O = numpy array of shape (len(row_label_unique),len(column_label_unique))
        row_labels_O = row labels of data_O
        column_labels_O = column labels of data_O

        OPTIONAL OUTPUT:
        row_variables_O = {"row_variables_I[0]:[...],..."} where each list is of len(row_labels_O)
        column_variables_O = {"row_variables_I[0]:[...],..."} where each list is of len(column_labels_O)
        '''
        pass;

    def extract_arrayFromListDict(self,data_I,row_label_I,column_label_I):
        '''convert a list of dictionary rows to a numpy array
        INPUT:
        
        
        '''
        pass;

    def convert_factorVector2responseMatrix(self,factors_I):
        '''convert a list of factors to a response matrix
        NOTES: required to convert pls to pls-da

        INPUT:
        factors_I = [] of strings
        
        OUTPUT:
        response_O = binary matrix of shape len(factor_I),len(factor_unique)
            where a 1 indicates association with the factor specified by that column
            and a 0 indicates no association with the factor specified by that column
        factors_O = list of unique factors of (the column names of response_O)
        '''
        pass;


