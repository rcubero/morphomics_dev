import morphomics
import numpy as np
import pandas as pd
import os
import umap
from sklearn.decomposition import PCA

class Protocols(object):
    
    def __init__(self, parameters, Parameters_ID):
        self.parameters = parameters
        self.file_prefix = "Morphomics.PID%d" % (
            int(Parameters_ID)
        )
        self.morphoframe = {}
        
        print(self.file_prefix)



    def Input(self):
        params = self.parameters["Input"]
        self.file_prefix = "%s.TMD-%s"%(self.file_prefix, params["barcode_filter"])
        
        # initialize output filename
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
            
        print(save_filename)

        # load the data
        self.morphoframe[params["morphoframe_name"]] = morphomics.io.load_data(
            folder_location=params["data_location_filepath"],
            extension=params["extension"],
            barcode_filter=params["barcode_filter"],
            save_filename=save_filename,
            conditions=params["conditions"],
            separated_by=params["separated_by"],
        )



    def Load_data(self):
        params = self.parameters["Load_data"]
        
        _morphoframe = {}
        for _c in params["conditions_to_include"]:
            print("...loading %s" % _c)
            save_filename = "%s/%s.%s-%s" % (params["folderpath_to_data"], params["filename_prefix"], params["separated_by"], _c)
            _morphoframe[_c] = morphomics.utils.load_obj(save_filename)

        self.morphoframe[params["morphoframe_name"]] = pd.concat([_morphoframe[_c] for _c in params["conditions_to_include"]], ignore_index=True)
        
        
        
    def Clean_frame(self):
        params = self.parameters["Clean_frame"]
        self.file_prefix = "%s.Cleaned"%(self.file_prefix)

        # initialize morphoframe to clean
        if params["morphoframe_filepath"]:
            _morphoframe = morphomics.utils.load_obj(params["morphoframe_filepath"])
        else:
            _morphoframe = self.morphoframe[params["morphoframe_name"]]

        _morphoframe = _morphoframe.loc[~_morphoframe.Barcodes.isna()].reset_index(
            drop=True
        )

        barlength_cutoff = float(params["barcode_size_cutoff"])
        _morphoframe["Barcode_length"] = _morphoframe.Barcodes.apply(lambda x: len(x))
        _morphoframe = _morphoframe.query(
            "Barcode_length >= @barlength_cutoff"
            ).reset_index(drop=True)

        if len(params["combine_conditions"]) > 0:
            for _cond, _before, _after in params["combine_conditions"]:
                _morphoframe.loc[_morphoframe[_cond].isin(_before), _cond] = _after

        if len(params["restrict_conditions"]) > 0:
            for _cond, _restricts, _action in params["restrict_conditions"]:
                assert _cond in _morphoframe.keys(), "%s not in morphoframe..."%_cond
                if _action == "drop":
                    for _restrictions in _restricts:
                        print("Filtering out %s from %s..."%(_restrictions, _cond))
                        assert _restrictions in _morphoframe[_cond].unique(), "%s not in the available condition..."%_restrictions
                        _morphoframe = _morphoframe.loc[~_morphoframe[_cond].str.contains(_restrictions)].reset_index(drop=True)
                elif _action == "keep":
                    _restrictions = "|".join(_restricts)
                    print("Keeping %s in %s..."%(_restrictions, _cond))
                    _morphoframe = _morphoframe.loc[_morphoframe[_cond].str.contains(_restrictions)].reset_index(drop=True)
                else:
                    print("Warning for ", _cond, _restricts)
                    print("Third column must be either 'drop' or 'keep'...")
                    print("Nothing will be done...")

        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
            morphomics.utils.save_obj(_morphoframe, save_filename)
            
        self.morphoframe[params["morphoframe_name"]] = _morphoframe

        
        
    def Bootstrap(self):
        params = self.parameters["Bootstrap"]
        self.file_prefix = "%s.Bootstrapped"%(self.file_prefix)
        
        if params["morphoframe_filepath"]:
            _morphoframe = morphomics.utils.load_obj(params["morphoframe_filepath"])
        else:
            _morphoframe = self.morphoframe[params["morphoframe_name"]]    
            
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
            
        if params["ratio"] == 0:
            ratio = None
        
        bootstrapped_frame = (
            morphomics.Analysis.bootstrapping.get_subsampled_population_from_infoframe(
                _morphoframe,
                condition_column=params["condition_column"],
                bootstrap_conditions=params["bootstrap_conditions"],
                bootstrap_resolution=params["bootstrap_resolution"],
                N_pop=params["N_pop"],
                N_samples=params["N_samples"],
                rand_seed=params["rand_seed"],
                ratio=ratio,
                save_filename=save_filename,
            )
        )
        
        self.bootstrap_info = (
            bootstrapped_frame[params["bootstrap_resolution"]]
            .reset_index(drop=True)
            .astype("category")
        )
        
        self.morphoframe[params["bootstrapframe_name"]] = bootstrapped_frame
        
        if params["save_data"]:
            morphomics.utils.save_obj(self.morphoframe[params["bootstrapframe_name"]], save_filename)
            morphomics.utils.save_obj(self.bootstrap_info, "%s-BootstrapInfo"%save_filename)
            
            
            
    def Persistence_Images(self):
        params = self.parameters["Persistence_Images"]
        self.file_prefix = "%s.PersistenceImages"%(self.file_prefix)
        
        if params["morphoframe_filepath"]:
            _morphoframe = morphomics.utils.load_obj(params["morphoframe_filepath"])
        else:
            _morphoframe = self.morphoframe[params["morphoframe_name"]]    
            
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
        
        # this is feature that will be developed in the future
        if params["barcode_weight"] == 0:
            barcode_weight = None
            
        self.PI_matrix = morphomics.Analysis.reduction.get_images_array_from_infoframe(
            _morphoframe,
            xlims=params["xlims"],
            ylims=params["ylims"],
            bw_method=params["bw_method"],
            norm_method=params["norm_method"],
            barcode_weight=barcode_weight,
            save_filename=save_filename,  # save the persistence images
        )
        
        
        
    def UMAP(self):
        params = self.parameters["UMAP"]
        self.file_prefix = "%s.UMAP"%(self.file_prefix)
        
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
            
        if params["PersistenceImages_filepath"]:
            self.PI_matrix = morphomics.utils.load_obj(params["PersistenceImages_filepath"])
            
        if params["filter_pixels"]:
            if params["filteredpixelindex_filepath"]:
                _tokeep = morphomics.utils.load_obj(params["filteredpixelindex_filepath"])
            else:
                _tokeep = np.where(
                    np.std(self.PI_matrix, axis=0) >= params["pixel_std_cutoff"]
                )[0]
                
            self.PI_matrix = np.array([np.array(self.PI_matrix[_i][_tokeep]) for _i in np.arange(len(self.PI_matrix))])

            if params["save_data"]:
                morphomics.utils.save_obj(_tokeep, "%s-FilteredIndex" % (save_filename))
                morphomics.utils.save_obj(self.PI_matrix, "%s-FilteredMatrix" % (save_filename))
                
        if params["run_PCA"]:
            F_PCA = PCA(n_components=params["n_PCs"])
            self.PI_matrix = F_PCA.fit_transform(self.PI_matrix)
        
            if params["save_data"]:
                morphomics.utils.save_obj(F_PCA, "%s-PCAfunction" % (save_filename) )
                morphomics.utils.save_obj(self.PI_matrix, "%s-PCAcoords" % (save_filename) )

        F_umap = umap.UMAP(
            n_neighbors=params["n_neighbors"],
            min_dist=params["min_dist"],
            spread=params["spread"],
            random_state=params["random_state"],
            n_components=params["n_components"],
            metric=params["metric"],
            densmap=bool(params["densmap"]),
        )

        self.X_umap = F_umap.fit_transform(self.PI_matrix)

        if params["save_data"]:
            morphomics.utils.save_obj(F_umap, "%s-UMAPfunction%dD" % (save_filename, params["n_components"]) )
            morphomics.utils.save_obj(self.X_umap, "%s-UMAPcoords%dD" % (save_filename, params["n_components"]) )
        
    
    
    def Palantir(self):
        params = self.parameters["Palantir"]
        self.file_prefix = "%s.Palantir"%(self.file_prefix)
        
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
            
        if params["X_umap_filepath"]:
            self.X_umap = morphomics.utils.load_obj(params["X_umap_filepath"])
        
        distances = morphomics.reduction.palantir_diffusion_maps(
            self.X_umap, 
            n_components=params["n_diffusion_components"], 
            knn=params["knn_diffusion"],
        )

        self.fdl = morphomics.reduction.force_directed_layout(
            distances["kernel"], 
            random_seed=params["fdl_random_seed"]
        )

        if params["save_data"]:
            morphomics.utils.save_obj(distances, "%s-PalantirDistances" % (save_filename) )
            morphomics.utils.save_obj(self.fdl, "%s-PalantirFDCoords" % (save_filename) )
        
        
        
    def Prepare_csv(self):
        params = self.parameters["Prepare_csv"]
        
        if params["declare_filepaths"] == 0:
            assert 'X_umap' in self.keys(), "Run UMAP first!"
            assert 'bootstrap_info' in self.keys(), "Bootstrap_info is not found"
        else:
            self.X_umap = morphomics.utils.load_obj( params["UMAP_filepath"] )
            self.bootstrap_info = morphomics.utils.load_obj( params["BootstrapInfo_filepath"] )
        
        _bootstrap_info = self.bootstrap_info.copy()
        for dims in range(self.X_umap.shape[1]):
            _bootstrap_info["UMAP_%d"%(dims+1)] = X_umap[:, dims]
        _bootstrap_info.to_csv( params["OutputCsv_filepath"] )
        

            
        
    def Mapping(self):
        params = self.parameters["Mapping"]
        self.file_prefix = "%s.Mapping"%(self.file_prefix)
            
        if params["save_data"]:
            if params["file_prefix"] == 0:
                params["file_prefix"] = self.file_prefix
            if params["save_folder"] == 0:
                params["save_folder"] = os.getcwd()
            save_filename = "%s/%s" % (params["save_folder"], params["file_prefix"])
        else:
            save_filename = None
            
        if os.path.isfile(params["F_umap_filepath"]):
            F_umap = morphomics.utils.load_obj(params["F_umap_filepath"])
        else:
            print("!!! IMPORTANT !!!")
            print("Please provide the filepath to the UMAP function that generated the morphological spectrum.")
            exit()
            

        if params["PersistenceImages_filepath"]:
            self.PI_matrix = morphomics.utils.load_obj(params["PersistenceImages_filepath"])
            
        if params["filter_pixels"]:
            if os.path.isfile(params["FilteredPixelIndex_filepath"]):
                _tokeep = morphomics.utils.load_obj(params["filteredpixelindex_filepath"])
            else:
                print("!!! IMPORTANT !!!")
                print("It is important that the indices filtered in the persistence image is consistent with that in the generated morphological spectrum.")
                exit()
                
            self.PI_matrix = np.array([np.array(self.PI_matrix[_i][_tokeep]) for _i in np.arange(len(self.PI_matrix))])

            if params["save_data"]:
                morphomics.utils.save_obj(self.PI_matrix, "%s-FilteredMatrix" % (save_filename))
                
        if params["run_PCA"]:
            if os.path.isfile(params["F_PCA_filepath"]):
                F_PCA = morphomics.utils.load_obj(params["F_PCA_filepath"])
            else:
                print("!!! IMPORTANT !!!")
                print("It is important that the PCA function is consistent with that in the generated morphological spectrum.")
                exit()    
    
            self.PI_matrix = F_PCA.transform(self.PI_matrix)
            
            if params["save_data"]:
                morphomics.utils.save_obj(self.PI_matrix, "%s-PCAcoords" % (save_filename))
        
        self.X_umap = F_umap.transform(self.PI_matrix)
        
        if params["save_data"]:
            morphomics.utils.save_obj(self.X_umap, "%s-UMAPcoords%dD" % (save_filename, params["n_components"]) )
        
        
        
    def Clear_morphoframe(self):
        self.morphoframe = None