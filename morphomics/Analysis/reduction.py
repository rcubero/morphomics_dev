"""
morphomics : dimensionality reduction tools

Author: Ryan Cubero
"""
import numpy as np
import concurrent.futures


from scipy.linalg import svd
from scipy.spatial.distance import cdist, squareform
from scipy.sparse import csr_matrix, find, issparse, coo_matrix
from scipy.sparse.linalg import eigs

from morphomics.Topology import analysis
from morphomics.utils import save_obj
from morphomics.utils import norm_methods, scipy_metric
from morphomics.default_parameters import defaults

try:
    from fa2 import ForceAtlas2
except ImportError as exc:
    raise ImportError(
        "ForceAltas2 is not installed. " + "Please install it by doing: pip install fa2"
    ) from exc
    
try:
    import umap
except ImportError as exc:
    raise ImportError(
        "UMAP is not installed. " + "Please install it by doing: pip install umap-learn"
    ) from exc


def _get_persistence_image_data_single(ar):
    """
    ar[0]:   persistence barcode
    ar[1,2]: x, y-lims
    ar[3]:   bw-method
    ar[4]:   normalization method (see morphomics.utils.norm_methods)
    ar[5]:   bar weights
    """
    if len(ar[0]) >= 0:
        res = analysis.get_persistence_image_data(
            ar[0], xlims=ar[1], ylims=ar[2], bw_method=ar[3], norm_method=ar[4], weights=ar[5]
        )
    else:
        res = []
    return res, ar[0]


def _get_pairwise_distance_from_persistence(
    imgs1, metric="l1", chunks=10, to_squareform=True
):
    """
    Returns mean and spread (standard deviation) of the point cloud of persistence images
    """
    N = len(imgs1)
    distances = np.zeros((N, N))

    # think about chunking this portion of the code
    splits = np.array_split(imgs1, chunks)
    _index = np.array_split(np.arange(N), chunks)
    splits_index = np.hstack(
        [
            0,
            [
                len(np.hstack([_index[i] for i in np.arange(j + 1)]))
                for j in np.arange(len(splits))
            ],
        ]
    )

    for i in np.arange(len(splits)):
        for j in np.arange(i, len(splits)):
            distances[
                splits_index[i] : splits_index[i + 1],
                splits_index[j] : splits_index[j + 1],
            ] = cdist(splits[i], splits[j], metric=scipy_metric[metric])

    # since there are non-zero lower diagonal elements
    distances[np.tril_indices(N)] = 0.0
    # symmetrize distance matrix
    distances = distances + distances.T

    if to_squareform == True:
        return squareform(distances)
    else:
        return distances


def get_images_array(
    p1, xlims=None, ylims=None, bw_method=None, norm_method="sum", 
):
    """
    Computes persistence images and returns an Nx10000 array,
    N is the number of barcodes in p1
    """
    # get the birth and death distance limits for the persistence images
    _xlims, _ylims = analysis.get_limits(p1)
    if xlims is None:
        xlims = _xlims
    if ylims is None:
        ylims = _ylims

    imgs1 = []
    p1_lims = []
    for p in p1:
        p1_lims.append([p, xlims, ylims, bw_method, norm_method])

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # calculate the persistence images in parallel
        for x, y in executor.map(
            _get_persistence_image_data_single, p1_lims, chunksize=200
        ):
            imgs1.append(x)

    N = len(imgs1)
    images = []
    for i in np.arange(N):
        images.append(imgs1[i].flatten() / norm_methods[norm_method](imgs1[i]))

    return images


def get_images_array_from_infoframe(
    _info_frame,
    xlims=None,
    ylims=None,
    bw_method=None,
    norm_method="sum",
    barcode_weight=None,
    save_filename=None,
):
    """Computes persistence images from info_frame for all specified conditions

    Args:
        info_frame (DataFrame): dataframe containing the barcodes that will be transformed into persistence images.
        xlims ([float,float], optional): left and right limits of the persistence image. Defaults to None.
        ylims ([float,float], optional): bottom and top limits of the persistence image. Defaults to None.
        bw_method (float, optional): spread of the Gaussian KDE. Defaults to None.
        norm_method (str, optional): method for normalizing persistence image. Can be either ["sum", "min", "max", "ave", "std", "l1" and "l2"]. Defaults to "sum".
        barcode_weight (list of float, optional): Defaults to None.

    Returns:
        DataFrame: info_frame with persistence image as a new column
    """
    assert (
        "Barcodes" in _info_frame.keys()
    ), "Missing `Barcodes` column in info_frame..."

    print("Calculating persistence images...")

    p1 = _info_frame["Barcodes"]
        
    # get the birth and death distance limits for the persistence images
    _xlims, _ylims = analysis.get_limits(p1)
    if xlims is None:
        xlims = _xlims
    if ylims is None:
        ylims = _ylims

    imgs1 = []
    p1_lims = []
    
    
    for _ind in np.arange(len(p1)):
        if barcode_weight is not None: 
            weights = barcode_weight[_ind]
        else:
            weights = None
        p1_lims.append([p1[_ind], xlims, ylims, bw_method, norm_method, weights])

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # parallelized calculation of the persistence images
        for x, y in executor.map(
            _get_persistence_image_data_single, p1_lims, chunksize=200
        ):
            imgs1.append(x)

    N = len(imgs1)
    images = []
    for i in np.arange(N):
        if len(imgs1[i]) > 0:
            images.append(imgs1[i].flatten() / norm_methods[norm_method](imgs1[i]))
        else:
            images.append(np.nan)

    if save_filename is not None:
        save_obj(np.array(images), save_filename)

    print("...done! \n")

    return np.array(images)


def get_distance_array(
    p1,
    xlims=None,
    ylims=None,
    norm_method="sum",
    metric="l1",
    chunks=10,
    to_squareform=True,
):
    """
    Computes and outputs array of pre-computed distances for heirarchical clustering
    """
    imgs1 = get_images_array(
        p1,
        xlims=xlims,
        ylims=ylims,
        norm_method=norm_method,
    )
    distances = _get_pairwise_distance_from_persistence(
        imgs1, metric=metric, chunks=chunks, to_squareform=to_squareform
    )

    if to_squareform == True:
        return squareform(distances)
    else:
        return distances


# learns coordinates of a force directed layout
def force_directed_layout(
    affinity_matrix, verbose=False, iterations=500, random_seed=10
):
    """ "
    Function to compute force directed layout from the affinity_matrix

    :param affinity_matrix: Sparse matrix representing affinities between cells
    :param cell_names: pandas Series object with cell names
    :param verbose: Verbosity for force directed layout computation
    :param iterations: Number of iterations used by ForceAtlas
    :return: Pandas data frame representing the force directed layout

    Code taken from: https://github.com/dpeerlab/Harmony/blob/master/src/harmony/plot.py

    New features on this script:
    - Added random_seed as an input to control output
    """

    np.random.seed(random_seed)
    init_coords = np.random.random((affinity_matrix.shape[0], 2))

    forceatlas2 = ForceAtlas2(
        # Behavior alternatives
        outboundAttractionDistribution=False,
        linLogMode=False,
        adjustSizes=False,
        edgeWeightInfluence=1.0,
        # Performance
        jitterTolerance=1.0,
        barnesHutOptimize=True,
        barnesHutTheta=1.2,
        multiThreaded=False,
        # Tuning
        scalingRatio=2.0,
        strongGravityMode=False,
        gravity=1.0,
        # Log
        verbose=verbose,
    )

    positions = forceatlas2.forceatlas2(
        affinity_matrix, pos=init_coords, iterations=iterations
    )
    positions = np.array(positions)

    return positions


def _get_sparse_matrix_from_indices_distances_umap(
    knn_indices, knn_dists, n_obs, n_neighbors
):
    """
    Code taken from: https://github.com/scverse/scanpy/blob/1fbbfcdbb53dda7f6ccab60098a9f9cd141a8025/scanpy/neighbors/__init__.py#L346
    """
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()


def palantir_diffusion_maps(
    X,
    metric="manhattan",
    n_components=10,
    knn=30,
    alpha=0,
    seed=None,
    metric_kwds=None,
):
    """Run Diffusion maps using the adaptive anisotropic kernel
    :param X: if sparse matrix, must be the adjacency matrix, otherwise, PCA projections of the data
    :param metric:
    :param n_components: Number of diffusion components
    :param knn: Number of nearest neighbors for graph construction
    :param alpha: Normalization parameter for the diffusion operator
    :param seed: Numpy random seed, randomized if None, set to an arbitrary integer for reproducibility
    :param metric_kwds:
    :return: Diffusion components, corresponding eigen values and the diffusion operator

    Code taken from: https://github.com/dpeerlab/Palantir/blob/79956db130b64ad3ac441a85887a51bba1e09a12/src/palantir/utils.py

    New features on this script:
    - Pulled the UMAP implementations for the nearest neighbors
    """

    # Determine the kernel
    N_data, _ = X.shape

    if not issparse(X):
        print("Determing nearest neighbor graph...")
        from umap.umap_ import nearest_neighbors

        if seed is not None:
            _seed = seed + 1024
        else:
            _seed = seed

        knn_indices, knn_distances, _ = nearest_neighbors(
            X,
            knn,
            random_state=_seed,
            metric=metric,
            metric_kwds=metric_kwds,
            angular=False,
            verbose=False,
        )

        kNN = _get_sparse_matrix_from_indices_distances_umap(
            knn_indices, knn_distances, N_data, knn
        )

        # Adaptive k
        print("Calculating affinites using adaptive Gaussian kernels...")
        adaptive_k = int(np.floor(knn / 3))
        adaptive_std = np.zeros(N_data)

        for i in np.arange(N_data):
            adaptive_std[i] = np.sort(kNN.data[kNN.indptr[i] : kNN.indptr[i + 1]])[
                adaptive_k - 1
            ]

        # Kernel
        x, y, dists = find(kNN)

        # X, y specific stds
        dists = dists / adaptive_std[x]
        W = csr_matrix((np.exp(-dists), (x, y)), shape=[N_data, N_data])

        # Diffusion components
        kernel = W + W.T
    else:
        kernel = X

    # Markov
    print("Calculating Markov chain...")
    D = np.ravel(kernel.sum(axis=1))

    if alpha > 0:
        # L_alpha
        D[D != 0] = D[D != 0] ** (-alpha)
        mat = csr_matrix((D, (range(N_data), range(N_data))), shape=[N_data, N_data])
        kernel = mat.dot(kernel).dot(mat)
        D = np.ravel(kernel.sum(axis=1))

    D[D != 0] = 1 / D[D != 0]
    T = csr_matrix((D, (range(N_data), range(N_data))), shape=[N_data, N_data]).dot(
        kernel
    )

    print("Implementing eigenvalue decomposition...")
    # Eigenvalue dcomposition
    np.random.seed(seed)
    v0 = np.random.rand(min(T.shape))
    D, V = eigs(T, n_components, tol=1e-4, maxiter=1000, v0=v0)
    D = np.real(D)
    V = np.real(V)
    inds = np.argsort(D)[::-1]
    D = D[inds]
    V = V[:, inds]

    # Normalize
    for i in range(V.shape[1]):
        V[:, i] = V[:, i] / np.linalg.norm(V[:, i])

    # Create result
    res = {"T": T, "EigenVectors": V, "EigenValues": D}
    res["kernel"] = kernel
    print("Done!")

    return res


def calculate_umap(
    barcodes,
    image_parameters=None,
    UMAP_parameters=None,
    save_folder=None,
    save_prefix=None,
):
    """
    Calculates the UMAP representation of the distance matrix X
    """
    # use default if image_parameters is not given
    complete_keys = ["xlims", "ylims", "norm_method", "metric", "chunks", "cutoff"]

    if image_parameters is None:
        image_parameters = {}

    for keys in complete_keys:
        try:
            image_parameters[keys]
        except:
            print("image_parameters: %s is not given. Reverting to default" % keys)
            image_parameters[keys] = defaults["image_parameters"][keys]

    # calculate distance matrix
    print("Calculating distance matrix...")
    X = get_distance_array(
        barcodes,
        xlims=image_parameters["xlims"],
        ylims=image_parameters["ylims"],
        norm_method=image_parameters["norm_method"],
        metric=image_parameters["metric"],
        chunks=image_parameters["chunks"],
        barcode_size_cutoff=image_parameters["cutoff"],
        to_squareform=True,
    )

    # use default if UMAP_parameters is not given
    complete_keys = ["N_dims", "n_neighbors", "min_dist", "spread", "random_state"]

    if UMAP_parameters is None:
        UMAP_parameters = {}

    for keys in complete_keys:
        try:
            UMAP_parameters[keys]
        except:
            print("UMAP_parameters: %s is not given. Reverting to default" % keys)
            UMAP_parameters[keys] = defaults["UMAP_parameters"][keys]

    print("Calculating singular value decomposition...")
    U, _, _ = svd(X)
    X_reduced = U.T[0 : UMAP_parameters["N_dims"]]

    print("Calculating UMAP representation...")
    X_umap = umap.UMAP(
        n_neighbors=UMAP_parameters["n_neighbors"],
        min_dist=UMAP_parameters["min_dist"],
        spread=UMAP_parameters["spread"],
        random_state=UMAP_parameters["random_state"],
    ).fit_transform(X_reduced.T)

    if save_folder is not None:
        if save_prefix is None:
            save_prefix = ""
        print("Saving results...")
        save_obj(X_umap, "%s/UMAP_%s" % (save_folder, save_prefix))

    print("Done!")
    return X_umap