"""
morphomics : plotting tools

Author: Ryan Cubero
"""
import numpy as np
import pickle as pkl

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, convex_hull_plot_2d

try:
    import ipyvolume as ipv
except ImportError as exc:
    raise ImportError(
        "ipyvolume is not installed. " + "Please install it by doing: pip install ipyvolume"
    ) from exc


FAD_conditions = ["3mpos", "6mpos"]
Ckp25_conditions = ["1w", "2w", "6w"]
dev_conditions = ["P7", "P15", "P22"]


def plot_convex_hulled(
    X_reduced,
    foreground_regions,
    background_regions,
    brain_conds,
    brain_labels,
    conditions,
    sizes,
    colors,
    pre_cond,
    savefile,
):
    xmax, xmin = np.amax(X_reduced[:, 0]), np.amin(X_reduced[:, 0])
    ymax, ymin = np.amax(X_reduced[:, 1]), np.amin(X_reduced[:, 1])

    fig, ax = plt.subplots(dpi=300)
    fig.set_size_inches(12, 10)

    for labs in background_regions:
        for conds in conditions:
            inds = np.where(
                (np.array(brain_conds) == conds) * (np.array(brain_labels) == labs)
            )[0]
            # ax.scatter(X_reduced[inds][:,0], X_reduced[inds][:,1], s=30,
            #           c='lightgrey', marker='o', lw=1.5, alpha=0.4, rasterized=True)

            distances = pdist(X_reduced[inds])
            distances = squareform(distances)
            graph = (distances <= 1.0).astype("int")
            graph = csr_matrix(graph)
            n_components, labels = connected_components(
                csgraph=graph, directed=False, return_labels=True
            )

            inds = inds[np.where(labels == 0)[0]]
            hull = ConvexHull(X_reduced[inds])
            for simplex in hull.simplices:
                plt.plot(
                    X_reduced[inds][simplex, 0],
                    X_reduced[inds][simplex, 1],
                    color=colors[conds],
                    ls="--",
                    lw=2,
                )

    for cond in conditions:
        for labs in foreground_regions:
            inds = np.where(
                (np.array(brain_conds) == cond) * (np.array(brain_labels) == labs)
            )[0]
            if pre_cond != None:
                label = "%s, %s %s" % (labs, pre_cond, cond)
            else:
                label = "%s, %s" % (labs, cond)
            ax.scatter(
                X_reduced[:, 0][inds],
                X_reduced[:, 1][inds],
                s=sizes[cond],
                c=colors[cond],
                marker="o",
                edgecolor="k",
                lw=0.4,
                label=label,
                rasterized=True,
            )

    ax.set_xlabel("UMAP 1", fontsize=14)
    ax.set_ylabel("UMAP 2", fontsize=14)

    ax.set_xlim(left=xmin * (1.1), right=xmax * (1.1))
    ax.set_ylim(bottom=ymin * (1.1), top=ymax * (1.1))

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    ax.legend(loc="lower right", fontsize=12)

    plt.savefig(savefile, bbox_inches="tight", dpi=300)


def plot_convex_hulled_MF_dev(
    X_reduced,
    foreground_region,
    brain_conds,
    brain_labels,
    brain_genders,
    conditions,
    sizes,
    colors_M,
    colors_F,
    pre_cond,
    savefile,
):
    xmax, xmin = np.amax(X_reduced[:, 0]), np.amin(X_reduced[:, 0])
    ymax, ymin = np.amax(X_reduced[:, 1]), np.amin(X_reduced[:, 1])

    marker = {}
    marker["M"] = "o"
    marker["F"] = "^"

    linestyle = {}
    linestyle["M"] = "--"
    linestyle["F"] = "dotted"

    size_offset = {}
    size_offset["M"] = 0
    size_offset["F"] = 20

    def plot_convex_hulled_sex_dev(foreground_gender, background_gender):
        fig, ax = plt.subplots(dpi=300)
        fig.set_size_inches(12, 10)

        for conds in conditions:
            inds = np.where(
                (np.array(brain_conds) == conds)
                * (np.array(brain_labels) == foreground_region)
                * (np.array(brain_genders) == background_gender)
            )[0]

            distances = pdist(X_reduced[inds])
            distances = squareform(distances)
            graph = (distances <= 0.8).astype("int")
            graph = csr_matrix(graph)
            n_components, labels = connected_components(
                csgraph=graph, directed=False, return_labels=True
            )

            largest_component = np.argmax(
                [len(np.where(labels == i)[0]) for i in np.unique(labels)]
            )
            inds = inds[np.where(labels == largest_component)[0]]
            hull = ConvexHull(X_reduced[inds])

            if background_gender == "M":
                color = colors_M[conds]
            elif background_gender == "F":
                color = colors_F[conds]
            else:
                color = "darkgrey"

            for simplex in hull.simplices:
                plt.plot(
                    X_reduced[inds][simplex, 0],
                    X_reduced[inds][simplex, 1],
                    color=color,
                    ls=linestyle[background_gender],
                    lw=2,
                )

        for conds in conditions:
            inds = np.where(
                (np.array(brain_conds) == conds)
                * (np.array(brain_labels) == foreground_region)
                * (np.array(brain_genders) == foreground_gender)
            )[0]
            if pre_cond != None:
                label = "%s, %s %s (%s)" % (
                    foreground_region,
                    pre_cond,
                    conds,
                    foreground_gender,
                )
            else:
                label = "%s, %s (%s)" % (foreground_region, conds, foreground_gender)

            if foreground_gender == "M":
                color = colors_M[conds]
            elif foreground_gender == "F":
                color = colors_F[conds]
            else:
                color = "darkgrey"

            ax.scatter(
                X_reduced[:, 0][inds],
                X_reduced[:, 1][inds],
                s=sizes[conds] + size_offset[foreground_gender],
                c=color,
                marker=marker[foreground_gender],
                edgecolor="k",
                lw=0.4,
                label=label,
                rasterized=True,
            )

        ax.set_xlabel("UMAP 1", fontsize=14)
        ax.set_ylabel("UMAP 2", fontsize=14)

        ax.set_xlim(left=xmin * (1.1), right=xmax * (1.1))
        ax.set_ylim(bottom=ymin * (1.1), top=ymax * (1.1))

        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        ax.legend(loc="lower right", fontsize=12)

        plt.savefig(
            savefile + "_%s.pdf" % foreground_gender, bbox_inches="tight", dpi=300
        )

    plot_convex_hulled_sex_dev("M", "F")
    plot_convex_hulled_sex_dev("F", "M")


def plot_convex_hulled_MF_deg_spatial(
    X_reduced,
    foreground_region,
    brain_conds,
    brain_labels,
    brain_genders,
    conditions,
    sizes,
    colors_M,
    colors_F,
    pre_cond,
    savefile,
):
    xmax, xmin = np.amax(X_reduced[:, 0]), np.amin(X_reduced[:, 0])
    ymax, ymin = np.amax(X_reduced[:, 1]), np.amin(X_reduced[:, 1])

    marker = {}
    marker["M"] = "o"
    marker["F"] = "^"

    linestyle = {}
    linestyle["M"] = "--"
    linestyle["F"] = "dotted"

    size_offset = {}
    size_offset["M"] = 0
    size_offset["F"] = 20

    def plot_convex_hulled_sex_deg(foreground_gender, background_gender):
        fig, ax = plt.subplots(dpi=300)
        fig.set_size_inches(12, 10)

        inds = np.where(np.array(brain_labels) != foreground_region)[0]
        ax.scatter(
            X_reduced[:, 0][inds],
            X_reduced[:, 1][inds],
            s=30,
            c="lightgrey",
            marker="o",
            alpha=0.4,
            rasterized=True,
        )

        for conds in conditions:
            inds = np.where(
                (np.array(brain_conds) == conds)
                * (np.array(brain_labels) == foreground_region)
                * (np.array(brain_genders) == background_gender)
            )[0]

            distances = pdist(X_reduced[inds])
            distances = squareform(distances)
            graph = (distances <= 0.8).astype("int")
            graph = csr_matrix(graph)
            n_components, labels = connected_components(
                csgraph=graph, directed=False, return_labels=True
            )

            largest_component = np.argmax(
                [len(np.where(labels == i)[0]) for i in np.unique(labels)]
            )
            inds = inds[np.where(labels == largest_component)[0]]
            hull = ConvexHull(X_reduced[inds])

            if background_gender == "M":
                color = colors_M[conds]
            elif background_gender == "F":
                color = colors_F[conds]
            else:
                color = "darkgrey"

            for simplex in hull.simplices:
                plt.plot(
                    X_reduced[inds][simplex, 0],
                    X_reduced[inds][simplex, 1],
                    color=color,
                    ls=linestyle[background_gender],
                    lw=2,
                )

        for conds in conditions:
            inds = np.where(
                (np.array(brain_conds) == conds)
                * (np.array(brain_labels) == foreground_region)
                * (np.array(brain_genders) == foreground_gender)
            )[0]

            if conds in FAD_conditions:
                label = "%s, 5xFAD %s (%s)" % (
                    foreground_region,
                    conds,
                    foreground_gender,
                )
            elif conds in Ckp25_conditions:
                label = "%s, Ckp25 %s (%s)" % (
                    foreground_region,
                    conds,
                    foreground_gender,
                )
            elif conds in dev_conditions:
                label = "%s, Dev %s (%s)" % (
                    foreground_region,
                    conds,
                    foreground_gender,
                )
            else:
                label = "%s, %s (%s)" % (foreground_region, conds, foreground_gender)

            if foreground_gender == "M":
                color = colors_M[conds]
            elif foreground_gender == "F":
                color = colors_F[conds]
            else:
                color = "darkgrey"

            ax.scatter(
                X_reduced[:, 0][inds],
                X_reduced[:, 1][inds],
                s=sizes[conds] + size_offset[foreground_gender],
                c=color,
                marker=marker[foreground_gender],
                edgecolor="k",
                lw=0.4,
                label=label,
                rasterized=True,
            )

        ax.set_xlabel("UMAP 1", fontsize=14)
        ax.set_ylabel("UMAP 2", fontsize=14)

        ax.set_xlim(left=xmin * (1.1), right=xmax * (1.1))
        ax.set_ylim(bottom=ymin * (1.1), top=ymax * (1.1))

        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

        ax.legend(loc="lower right", fontsize=12)

        if savefile != "FALSE":
            plt.savefig(
                savefile + "_%s.pdf" % foreground_gender, bbox_inches="tight", dpi=300
            )

    plot_convex_hulled_sex_deg("M", "F")
    plot_convex_hulled_sex_deg("F", "M")


def quick_pyvol_UMAP(_X_umap, _info_frame, _conditions, _regions, _color):
    ipv.figure()
    ipv.scatter(
        _X_umap[:, 0],
        _X_umap[:, 1],
        _X_umap[:, 2],
        size=0.1,
        color="lightgrey",
        marker="sphere",
    )
    for region in _regions:
        _inds = np.where(
            (_info_frame.Region == region)
            * (_info_frame.Model.str.contains(_conditions))
        )[0]
        if len(_inds) > 0:
            ipv.scatter(
                _X_umap[:, 0][_inds],
                _X_umap[:, 1][_inds],
                _X_umap[:, 2][_inds],
                size=0.5,
                color=_color[region],
                marker="sphere",
                label=region,
            )
    ipv.show()
    

def _get_condition(_model):
    if _model in ["Cx3cr1_het", "Light_rearing"]:
        return "Development"
    if _model in ["rd1", "rd10"]:
        return "Degeneration"


def _pyvol_scatter_UMAP(_X_umap, _info_frame, _model, _region, _cmaps):
    ipv.scatter(
        _X_umap[:, 0],
        _X_umap[:, 1],
        _X_umap[:, 2],
        size=1,
        color="lightgrey",
        marker="sphere",
        alpha=0.3,
    )

    query_frame = _info_frame.query("Region == @_region & Model == @_model")
    timepoints = query_frame.Time.unique()
    time_vals = [
        int(_time.split("_")[0][1:]) if _time.split("_")[0][0] == "P" else 30
        for _time in timepoints
    ]
    _cond = _get_condition(_model)

    for _t in np.arange(len(timepoints)):
        _time = timepoints[_t]
        _inds = query_frame.query("Time == @_time").index.values

        color = [
            _cmaps[_cond][_region]((time_vals[_t] + 5) / (max(time_vals) * 1.01))
        ] * len(_inds)

        if len(_inds) > 0:
            ipv.scatter(
                _X_umap[:, 0][_inds],
                _X_umap[:, 1][_inds],
                _X_umap[:, 2][_inds],
                size=2,
                color=color,
                marker="sphere",
                label=_time,
                alpha=0.7,
            )


def quick_pyvol_UMAP_trajectory(_X_umap, _info_frame, _model, _regions, _cmaps):
    ipv.figure()
    ipv.style.box_off()
    for _region in _regions:
        _pyvol_scatter_UMAP(_X_umap, _info_frame, _model, _region, _cmaps)
    ipv.show()