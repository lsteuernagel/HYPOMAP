# from: https://github.com/lhqing/ALLCools/blob/master/ALLCools/clustering/ConsensusClustering.py 
# on 19.01.24 , last change with commit: 7c16590

import warnings
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

import anndata
import joblib
import leidenalg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from imblearn.ensemble import BalancedRandomForestClassifier
from natsort import natsorted
from scanpy.neighbors import Neighbors
from sklearn.metrics import (
    adjusted_rand_score,
    balanced_accuracy_score,
    confusion_matrix,
    pairwise_distances,
)
from sklearn.model_selection import cross_val_predict

# from ..plot import categorical_scatter


def _adata_to_coord_data(adata, coord_base):
    coord_data = pd.DataFrame(
        {f"{coord_base}_0": adata.obsm[f"X_{coord_base}"][:, 0], f"{coord_base}_1": adata.obsm[f"X_{coord_base}"][:, 1]}
    )
    return coord_data


def _r1_normalize(cmat):
    """
    Perofrm R1 normalization on the confusion matrix.

    Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

    Normalize the confusion matrix based on the total number of cells in each class
    x(i,j) = max(cmat(i,j)/diagnol(i),cmat(j,i)/diagnol(j))
    confusion rate between i and j is defined by the maximum ratio i is confused as j or j is confused as i.

    Parameters
    ----------
    cmat :
        the confusion matrix

    Returns
    -------
    the R1 normalized confusion matrix
    """
    dmat = cmat
    smat = np.diag(dmat) + 1  # in case some label has no correct prediction (0 in diag)
    dim = cmat.shape[0]
    xmat = np.zeros([dim, dim])
    for i in range(dim):
        for j in range(i + 1, dim):
            xmat[i, j] = xmat[j, i] = max(dmat[i, j] / smat[j], dmat[j, i] / smat[i])

    # scale matrix to 0-1
    xmat = xmat / np.max(xmat)

    return xmat


def _r2_normalize(cmat):
    """
    Perofrm R2 normalization on the confusion matrix.

    Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

    Normalize the confusion matrix based on the total number of cells.
    x(i,j) = max(cmat(i,j)+cmat(j,i)/N)
    N is total number of cells analyzed.
    Confusion rate between i and j is defined by the sum of i confused as j or j confused as i.
    Then divide by total number of cells.

    Parameters
    ----------
    cmat :
        the confusion matrix

    Returns
    -------
    the R2 normalized confusion matrix
    """
    emat = np.copy(cmat)
    s = np.sum(cmat)
    emat = emat + emat.T
    np.fill_diagonal(emat, 0)
    emat = emat * 1.0 / s

    # scale matrix to 0-1
    emat = emat / np.max(emat)

    return emat


def _leiden_runner(g, random_states, partition_type, **partition_kwargs):
    """
    Run leiden on a graph and return the partition.

    The leiden clustering repeated len(random_states) times with different random states,
    return all clusters as a pd.DataFrame.
    """
    results = []
    for seed in random_states:
        part = leidenalg.find_partition(g, partition_type, seed=seed, **partition_kwargs)
        groups = np.array(part.membership)
        groups = pd.Categorical(
            values=groups.astype("U"),
            categories=natsorted(np.unique(groups).astype("U")),
        )
        results.append(groups)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        result_df = pd.DataFrame(results, columns=random_states)
    return result_df


def _split_train_test_per_group(x, y, frac, max_train, random_state):
    """Split train test for each cluster and make sure there are enough cells for train."""
    y_series = pd.Series(y)
    # split train test per group
    train_idx = []
    test_idx = []
    outlier_idx = []
    for cluster, sub_series in y_series.groupby(y_series):
        if (cluster == -1) or (sub_series.size < 3):
            outlier_idx += sub_series.index.tolist()
        else:
            n_train = max(1, min(max_train, int(sub_series.size * frac)))
            is_train = sub_series.index.isin(sub_series.sample(n_train, random_state=random_state).index)
            train_idx += sub_series.index[is_train].tolist()
            test_idx += sub_series.index[~is_train].tolist()
    x_train = x[train_idx]
    y_train = y[train_idx]
    x_test = x[test_idx]
    y_test = y[test_idx]
    return x_train, y_train, x_test, y_test


def single_supervise_evaluation(clf, x_train, y_train, x_test, y_test, r1_norm_step=0.05, r2_norm_step=0.05):
    """Run supervise evaluation on confusion matrix."""
    # fit model
    clf.fit(x_train, y_train)

    # calc accuracy
    y_train_pred = clf.predict(x_train)
    accuracy_train = balanced_accuracy_score(y_true=y_train, y_pred=y_train_pred)
    print(f"Balanced accuracy on the training set: {accuracy_train:.3f}")
    y_test_pred = clf.predict(x_test)
    accuracy_test = balanced_accuracy_score(y_true=y_test, y_pred=y_test_pred)
    print(f"Balanced accuracy on the hold-out set: {accuracy_test:.3f}")

    # get confusion matrix
    y_pred = clf.predict(x_test)
    cmat = confusion_matrix(y_test, y_pred)

    # normalize confusion matrix
    r1_cmat = _r1_normalize(cmat)
    r2_cmat = _r2_normalize(cmat)
    m1 = np.max(r1_cmat)
    if np.isnan(m1):
        m1 = 1.0
    m2 = np.max(r2_cmat)

    cluster_map = {}
    while (len(cluster_map) == 0) and (m1 > 0) and (m2 > 0):
        m1 -= r1_norm_step
        m2 -= r2_norm_step

        # final binary matrix to calculate which clusters need to be merged
        judge = np.maximum.reduce([(r1_cmat > m1), (r2_cmat > m2)])
        if judge.sum() > 0:
            rows, cols = np.where(judge)
            edges = zip(rows.tolist(), cols.tolist())
            g = nx.Graph()
            g.add_edges_from(edges)
            for comp in nx.connected_components(g):
                to_label = comp.pop()
                for remain in comp:
                    cluster_map[remain] = to_label
    return clf, accuracy_test, cluster_map, cmat, r1_cmat, r2_cmat


class ConsensusClustering:
    def __init__(
        self,
        model=None,
        n_neighbors=25,
        metric="euclidean",
        min_cluster_size=10,
        leiden_repeats=200,
        leiden_resolution=1,
        target_accuracy=0.95,
        consensus_rate=0.7,
        random_state=0,
        train_frac=0.5,
        train_max_n=500,
        max_iter=50,
        n_jobs=-1,
    ):
        """
        Perform consensus clustering by multi-leiden clustering + supervised model evaluation.

        Parameters
        ----------
        model
            A supervised ML model, if not provided, will use the BalancedRandomForestClassifier from imblearn
        n_neighbors
            K of the KNN graph
        metric
            metric of the KNN graph
        min_cluster_size
            Minimum final cluster size to report
        consensus_rate
            Cutoff for the initial cluster separation from multi-leiden run
        leiden_repeats
            Repeat leiden clustering with different random states this number of times
        leiden_resolution
            Resolution parameter of leiden clustering
        random_state
            Overall random state to assure reproducibility
        train_frac
            fraction of cells per cluster used in training supervised model
        train_max_n
            maximum number of cells per cluster used in training supervised model, this is to prevent some
            large cluster dominant the training dataset. The actual cells used in training is the smaller one
            of `train_max_n` and `train_frac * cluster_size`
        max_iter
            maximum iteration of the cluster merge process
        n_jobs
            number of cpus
        """
        # input metrics
        self.min_cluster_size = min_cluster_size
        self.consensus_rate = consensus_rate  # this prevents merging gradient clusters
        self.leiden_repeats = leiden_repeats
        self.leiden_resolution = leiden_resolution
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.n_neighbors = n_neighbors
        self.knn_metric = metric
        self.train_frac = train_frac
        self.train_max_n = train_max_n
        self.max_iter = max_iter
        self.n_obs, self.n_pcs = None, None
        self.X = None
        self._neighbors = None
        self.step_data = OrderedDict()
        self.target_accuracy = target_accuracy

        # multiple leiden clustering
        self.leiden_result_df = None
        self._multi_leiden_clusters = None

        # model training and outlier rescue
        self.supervise_model = model
        self._label_with_leiden_outliers = None
        self.label = None
        self.label_proba = None
        self.cv_predicted_label = None
        self.final_accuracy = None
        return

    def add_data(self, x):
        # just add X, but not doing computation
        # use this for step-by-step computation
        self.n_obs, self.n_pcs = x.shape
        self.X = x

    def fit_predict(self, x, leiden_kwds=None):
        self.add_data(x)

        # Construct KNN graph
        print("Computing nearest neighbor graph")
        self.compute_neighbors()

        # repeat Leiden clustering with different random seeds
        # summarize the results and determine initial clusters by hamming distance between leiden runs
        print("Computing multiple clustering with different random seeds")
        kwds = {}
        if leiden_kwds is None:
            leiden_kwds = kwds
        else:
            leiden_kwds.update(kwds)
        self.multi_leiden_clustering(**leiden_kwds)

        # merge the over clustering version by supervised learning
        self.supervise_learning()

        # assign outliers
        self.final_evaluation()

    def compute_neighbors(self):
        """Calculate KNN graph"""
        # nearest neighbors graph
        adata = anndata.AnnData(
            X=None,
            obs=pd.DataFrame([], index=[f"obs{i}" for i in range(self.n_obs)]),
            var=pd.DataFrame([], index=[f"var{i}" for i in range(self.n_pcs)]),
        )
        adata.obsm["X_pca"] = self.X
        # here neighbors should only use PCs
        self._neighbors = Neighbors(adata=adata)
        self._neighbors.compute_neighbors(
            n_neighbors=self.n_neighbors,
            knn=True,
            n_pcs=self.n_pcs,
            use_rep="X_pca",
            method="umap",
            metric=self.knn_metric,
            random_state=self.random_state,
        )
        return

    def multi_leiden_clustering(
        self,
        partition_type=None,
        partition_kwargs=None,
        use_weights=True,
        n_iterations=-1,
    ):
        """Run multiple leiden clustering with different random seeds and summarize the results."""
        if self._neighbors is None:
            raise ValueError("Run compute_neighbors first before multi_leiden_clustering")

        # convert neighbors to igraph
        g = self._neighbors.to_igraph()

        # generate n different seeds for each single leiden partition
        np.random.seed(self.random_state)
        leiden_repeats = self.leiden_repeats
        n_jobs = self.n_jobs
        random_states = np.random.choice(range(99999), size=leiden_repeats, replace=False)
        step = max(int(leiden_repeats / n_jobs), 10)
        random_state_chunks = [random_states[i : min(i + step, leiden_repeats)] for i in range(0, leiden_repeats, step)]

        results = []
        print(f"Repeating leiden clustering {leiden_repeats} times")
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            future_dict = {}
            for random_state_chunk in random_state_chunks:
                # flip to the default partition type if not over writen by the user
                if partition_type is None:
                    partition_type = leidenalg.RBConfigurationVertexPartition
                # prepare find_partition arguments as a dictionary, appending to whatever the user provided
                # it needs to be this way as this allows for the accounting of a None resolution
                # (in the case of a partition variant that doesn't take it on input)
                if partition_kwargs is None:
                    partition_kwargs = {}
                else:
                    if "seed" in partition_kwargs:
                        print("Warning: seed in the partition_kwargs will be ignored, use seed instead.")
                        del partition_kwargs["seed"]
                if use_weights:
                    partition_kwargs["weights"] = np.array(g.es["weight"]).astype(np.float64)
                partition_kwargs["n_iterations"] = n_iterations
                partition_kwargs["resolution_parameter"] = self.leiden_resolution
                # clustering proper
                future = executor.submit(
                    _leiden_runner,
                    g=g,
                    random_states=random_state_chunk,
                    partition_type=partition_type,
                    **partition_kwargs,
                )
                future_dict[future] = random_state_chunks

            for future in as_completed(future_dict):
                _ = future_dict[future]
                try:
                    data = future.result()
                    results.append(data)
                except Exception as exc:
                    print(f"_leiden_runner generated an exception: {exc}")
                    raise exc
        total_result = pd.concat(results, axis=1, sort=True)
        self.leiden_result_df = total_result
        cluster_count = self.leiden_result_df.apply(lambda i: i.unique().size)
        print(
            f"Found {cluster_count.min()} - {cluster_count.max()} clusters, "
            f"mean {cluster_count.mean():.1f}, std {cluster_count.std():.2f}"
        )
        # create a over-clustering version based on all the leiden runs
        print("Summarizing multiple clustering results")
        self._summarize_multi_leiden()
        return

    def _summarize_multi_leiden(self):
        """
        Summarize the multi_leiden results.

        Generate a raw cluster version simply based on the hamming distance
        between cells and split cluster with cutoff (consensus_rate)
        """
        # data: row is leiden run, column is cell
        data = self.leiden_result_df.T

        # group cell into raw clusters if their hamming distance < 1 - consensus_rate
        cur_cluster_id = 0
        clusters = {}
        while data.shape[1] > 1:
            seed_cell = data.pop(data.columns[0])
            distance = pairwise_distances(X=data.T, Y=seed_cell.values[None, :], metric="hamming").ravel()

            judge = distance < (1 - self.consensus_rate)
            this_cluster_cells = [seed_cell.name] + data.columns[judge].to_list()
            for cell in this_cluster_cells:
                clusters[cell] = cur_cluster_id
            data = data.loc[:, ~judge].copy()
            cur_cluster_id += 1
        if data.shape[1] == 1:
            # if there is only one cell remain
            clusters[data.columns[0]] = cur_cluster_id
        clusters = pd.Series(clusters).sort_index()

        # renumber clusters based on cluster size
        counts = clusters.value_counts()
        cluster_map = {c: i for i, c in enumerate(counts.index)}
        clusters = clusters.map(cluster_map)
        # renumber small clusters as -1
        counts = clusters.value_counts()
        small_clusters = counts[counts < self.min_cluster_size].index
        clusters[clusters.isin(small_clusters)] = -1

        print(f"{(clusters != -1).sum()} cells assigned to {clusters.unique().size - 1} raw clusters")
        print(f"{(clusters == -1).sum()} cells are multi-leiden outliers")
        self._multi_leiden_clusters = clusters.values
        return

    def _create_model(self, n_estimators=1000):
        """Init default model"""
        clf = BalancedRandomForestClassifier(
            n_estimators=n_estimators,
            criterion="gini",
            max_depth=None,
            min_samples_split=2,
            min_samples_leaf=2,
            min_weight_fraction_leaf=0.0,
            max_features="sqrt",
            max_leaf_nodes=None,
            min_impurity_decrease=0.0,
            bootstrap=True,
            oob_score=False,
            sampling_strategy="auto",
            replacement=False,
            n_jobs=self.n_jobs,
            random_state=self.random_state,
            verbose=0,
            warm_start=False,
            class_weight=None,
        )
        return clf

    def supervise_learning(self):
        """Perform supervised learning and cluster merge process"""
        if self._multi_leiden_clusters is None:
            raise ValueError(
                "Run multi_leiden_clustering first to get a " "clustering assignment before run supervise_learning."
            )

        n_cluster = np.unique(self._multi_leiden_clusters[self._multi_leiden_clusters != -1]).size
        if n_cluster == 1:
            print("There is only one cluster except for outliers, can not train supervise model on that.")
            self.label = np.zeros(self.n_obs, dtype=int)
            return
        print("\n=== Start supervise model training and cluster merging ===")

        x = self.X
        cur_y = self._multi_leiden_clusters.copy()
        score = None
        step = 0.1

        if self.supervise_model is None:
            # create default model if no model provided
            clf = self._create_model(n_estimators=500)
        else:
            clf = self.supervise_model
        for cur_iter in range(1, self.max_iter + 1):
            print(f"\n=== iteration {cur_iter} ===")
            n_labels = np.unique(cur_y[cur_y != -1]).size
            print(f"{n_labels} non-outlier labels")
            if n_labels < 2:
                print(f"Stop iteration because only {n_labels} cluster remain.")
                break

            x_train, y_train, x_test, y_test = _split_train_test_per_group(
                x=x,
                y=cur_y,
                frac=self.train_frac,
                max_train=self.train_max_n,
                random_state=self.random_state + cur_iter,
                # every time train-test split got a different random state
            )
            (
                clf,
                score,
                cluster_map,
                cmat,
                r1_cmat,
                r2_cmat,
            ) = single_supervise_evaluation(
                clf,
                x_train,
                y_train,
                x_test,
                y_test,
                r1_norm_step=step,
                r2_norm_step=step,
            )
            step = min(0.2, max(0.05, 2 * (self.target_accuracy - score)))

            # save step data for plotting
            self.step_data[cur_iter] = [
                cur_y,
                cmat,
                r1_cmat,
                r2_cmat,
                cluster_map,
                score,
                step,
            ]

            if score > self.target_accuracy:
                print(
                    f"Stop iteration because current accuracy {score:.3f}"
                    f" > target accuracy {self.target_accuracy:.3f}."
                )
                break

            # judge results
            if len(cluster_map) > 0:
                print(f"Merging {len(cluster_map)} clusters.")
                cur_y = pd.Series(cur_y).apply(lambda i: cluster_map[i] if i in cluster_map else i)
                # renumber labels from large to small
                ordered_map = {c: i for i, c in enumerate(cur_y[cur_y != -1].value_counts().index)}
                cur_y = pd.Series(cur_y).apply(lambda i: ordered_map[i] if i in ordered_map else i).values
            else:
                print("Stop iteration because there is no cluster to merge")
                break
        else:
            print("Stop iteration because reaching maximum iteration.")
        self._label_with_leiden_outliers = cur_y
        self.label = cur_y
        self.supervise_model = clf
        self.final_accuracy = score
        return

    def final_evaluation(self):
        """Evaluate the final model"""
        print("\n=== Assign final labels ===")

        # skip if there is only one cluster
        n_cluster = len(set(self.label[self.label != -1]))
        if n_cluster < 2:
            print(f"Skip final evaluation because only {n_cluster} cluster label exist.")
            # name all cluster as c0
            self.label = np.zeros(self.label.size, dtype=int)
            self.cv_predicted_label = [f"c{label}" for label in self.label]
            self.label_proba = np.ones(self.label.size, dtype=int)
            self.final_accuracy = 1
        else:
            # predict outliers
            outlier_x = self.X[self.label == -1]
            outlier_idx = np.where(self.label == -1)[0]
            if len(outlier_idx) != 0:
                outlier_predict = pd.Series(self.supervise_model.predict(outlier_x), index=outlier_idx)
                for cell, pred_label in outlier_predict.items():
                    self.label[cell] = pred_label
            print(
                "Assigned all the multi-leiden clustering outliers into clusters "
                "using the prediction model from final clustering version."
            )

            # final evaluation of non-outliers using cross val predict
            final_predict_proba = cross_val_predict(
                self.supervise_model,
                self.X,
                y=self.label,
                method="predict_proba",
                n_jobs=self.n_jobs,
                verbose=0,
                cv=10,
            )
            final_predict = pd.Series(np.argmax(final_predict_proba, axis=1))
            final_cell_proba = pd.Series(np.max(final_predict_proba, axis=1))
            final_acc = balanced_accuracy_score(self.label, final_predict.values)
            print(f"Final ten-fold CV Accuracy on all the cells: {final_acc:.3f}")
            self.cv_predicted_label = [f"c{label}" for label in final_predict]
            self.label_proba = final_cell_proba.values
            self.final_accuracy = final_acc

        self.label = [f"c{label}" for label in self.label]
        return

    def save(self, output_path):
        """Save the model"""
        joblib.dump(self, output_path)

    # def plot_leiden_cases(self, coord_data, coord_base="umap", plot_size=3, dpi=300, plot_n_cases=4, s=3):
    #     """Show some leiden runs with the biggest different as measured by ARI"""
    #     if isinstance(coord_data, anndata.AnnData):
    #         coord_data = _adata_to_coord_data(coord_data, coord_base)
    # 
    #     # choose some most different leiden runs by rand index
    #     sample_cells = min(1000, self.leiden_result_df.shape[0])
    #     sample_runs = min(50, self.leiden_result_df.shape[1])
    #     use_df = self.leiden_result_df.sample(sample_cells, replace=False).T.sample(sample_runs)
    #     rand_index_rank = (
    #         pd.DataFrame(
    #             pairwise_distances(use_df, metric=adjusted_rand_score),
    #             index=use_df.index,
    #             columns=use_df.index,
    #         )
    #         .unstack()
    #         .sort_values()
    #     )
    #     plot_cases = set()
    #     for pairs in rand_index_rank[:10].index:
    #         plot_cases.add(pairs[0])
    #         plot_cases.add(pairs[1])
    #         if len(plot_cases) > plot_n_cases:
    #             break
    #     plot_cases = list(plot_cases)[:plot_n_cases]
    # 
    #     # plot
    #     plot_data = coord_data.copy()
    #     fig, axes = plt.subplots(figsize=(plot_n_cases * plot_size, plot_size), ncols=plot_n_cases, dpi=dpi)
    # 
    #     for case, ax in zip(plot_cases, axes):
    #         plot_data[f"Leiden {case}"] = self.leiden_result_df[case].values
    #         categorical_scatter(ax=ax, data=plot_data, coord_base=coord_base, hue=f"Leiden {case}", s=s)
    #     return fig, axes

    def plot_before_after(self, coord_data, coord_base="umap", plot_size=3, dpi=300):
        """Plot the raw clusters from multi-leiden and final clusters after merge"""
        if isinstance(coord_data, anndata.AnnData):
            coord_data = _adata_to_coord_data(coord_data, coord_base)
        if len(self.step_data) == 0:
            print("No merge step to plot")
            return
        plot_data = coord_data.copy()
        # initial clusters
        cur_y, *_ = self.step_data[1]
        fig, axes = plt.subplots(figsize=(2 * plot_size, plot_size), ncols=2, dpi=dpi)
        ax = axes[0]
        plot_data["cur_y"] = cur_y
        _ = categorical_scatter(
            data=plot_data,
            ax=ax,
            hue="cur_y",
            coord_base=coord_base,
            palette="tab20",
            text_anno="cur_y",
            show_legend=False,
        )
        ax.set(title="Initial Labels From\nMulti-leiden Clustering")
        ax = axes[1]
        plot_data["cur_y"] = self.label
        _ = categorical_scatter(
            data=plot_data,
            ax=ax,
            hue="cur_y",
            coord_base=coord_base,
            palette="tab20",
            text_anno="cur_y",
            show_legend=False,
        )
        ax.set(title="Final Labels After Merge")
        return

    def plot_steps(self, coord_data, coord_base="umap", plot_size=3, dpi=300):
        """Plot the supervised learning and merge steps"""
        if len(self.step_data) == 0:
            print("No merge step to plot")
            return
        if isinstance(coord_data, anndata.AnnData):
            plot_data = _adata_to_coord_data(coord_data, coord_base)
        else:
            plot_data = coord_data.copy()
        self.plot_before_after(coord_data, coord_base=coord_base, plot_size=plot_size, dpi=dpi)

        # initial clusters
        for i, step in enumerate(sorted(self.step_data.keys())):
            (
                cur_y,
                cmat,
                r1_cmat,
                r2_cmat,
                cluster_map,
                score,
                step_size,
            ) = self.step_data[step]
            fig, axes = plt.subplots(figsize=(4 * plot_size, plot_size), ncols=4, dpi=dpi)
            ax = axes[0]
            sns.heatmap(ax=ax, data=cmat, cbar=None, cmap="Reds")
            ax.set(title="Confusion Matrix", ylabel=f"Step {step}")
            ax = axes[1]
            sns.heatmap(ax=ax, data=r1_cmat, cbar=None, cmap="Reds")
            ax.set(title="R1 Norm.")
            ax = axes[2]
            sns.heatmap(ax=ax, data=r2_cmat, cbar=None, cmap="Reds")
            ax.set(title="R2 Norm.")
            ax = axes[3]
            if len(cluster_map) > 0:
                involved_clusters = set(list(cluster_map.keys()) + list(cluster_map.values()))
                cur_y = pd.Series(cur_y, index=plot_data.index)
                cur_y = cur_y.apply(lambda i: cluster_map[i] if i in cluster_map else i)
                # if not involved, mark as -1
                cur_y = cur_y.apply(lambda i: i if i in involved_clusters else -1)
                plot_data["cur_y"] = cur_y
                n_color = cur_y.unique().size
                colors = list(sns.color_palette("tab10", n_color))
                cmap = {c: colors.pop() if c != -1 else (0.9, 0.9, 0.9) for c in cur_y.unique()}
                categorical_scatter(
                    ax=ax,
                    data=plot_data,
                    coord_base=coord_base,
                    hue="cur_y",
                    palette=cmap,
                    s=3,
                )
                ax.set(title=f"Step {step} Cells After Merge")
            else:
                ax.axis("off")
        return

    def plot_merge_process(self, plot_size=3):
        """Plot the change of accuracy during merge"""
        if len(self.step_data) == 0:
            print("No merge step to plot")
            return
        plot_data = []
        for step, (cur_y, *_, score, _) in self.step_data.items():
            cur_n_cluster = len(set(cur_y)) - 1
            plot_data.append({"Step": step, "# of clusters": cur_n_cluster, "Score": score})
        plot_data = pd.DataFrame(plot_data)

        fig, axes = plt.subplots(
            figsize=(plot_size * 2, plot_size),
            ncols=2,
            dpi=300,
            constrained_layout=True,
        )
        ax = axes[0]
        sns.lineplot(ax=ax, data=plot_data, x="Step", y="# of clusters", color="gray")
        ax = axes[1]
        sns.lineplot(ax=ax, data=plot_data, x="Step", y="Score", color="gray")
        ax.axhline(y=self.target_accuracy, color="red", linewidth=0.8, linestyle="--")
        ax.set(ylim=(plot_data["Score"].min() - 0.03, min(1, self.target_accuracy + 0.01)))
        return


def select_confusion_pairs(true_label, predicted_label, ratio_cutoff=0.001):
    """
    Select cluster pairs that are confusing (ratio_cutoff) between true and predicted labels

    Parameters
    ----------
    true_label : true cell labels
    predicted_label : predicted cell labels
    ratio_cutoff : ratio of clusters cutoff to define confusion

    Returns
    -------
    confused_pairs :
        list of cluster pair tuples
    """
    labels = pd.DataFrame({"true": true_label, "pred": predicted_label})
    confusion_matrix = labels.groupby("true")["pred"].value_counts().unstack().fillna(0).astype(int)

    row_sum = confusion_matrix.sum(axis=1)
    row_norm = (confusion_matrix / row_sum.values[:, None]).unstack()
    row_pairs = row_norm[row_norm > ratio_cutoff].reset_index().iloc[:, :2]
    col_sum = confusion_matrix.sum(axis=0)
    col_norm = (confusion_matrix / col_sum.values[None, :]).unstack()
    col_pairs = col_norm[col_norm > ratio_cutoff].reset_index().iloc[:, :2]

    include_pairs = set()
    for _, s in pd.concat([row_pairs, col_pairs]).iterrows():
        a, b = s.sort_values()
        if a == b:
            continue
        include_pairs.add((a, b))
    return list(include_pairs)
