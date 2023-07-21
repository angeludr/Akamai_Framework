event_table = {
    # Opens file
    "open_fits_file":                 ("OpenFits",
                                       "opening_fits_file",
                                       "branch"),
    
    "branch":                         ("Branch",
                                       "getting_extensions",
                                        None),
    
    # Gets info for each extension, computes bias, and subtracts overscan region
    "get_header":                     ("GetHeader",
                                       "getting_fits_header",
                                       "get_binning"),
    
    "get_binning":                    ("GetBinning",
                                       "getting_binning_dimensions",
                                       "get_overscan"),
    
    "get_overscan":                   ("GetOverscan",
                                       "getting_overscan_regions",
                                       "get_bias"),
    
    "get_bias":                       ("GetBias",
                                       "getting_bias",
                                       "plot_bias"),
    
    "plot_bias":                      ("PlotBias",
                                       "plotting_bias",
                                       "bias_subtraction"),
    
    "bias_subtraction":               ("BiasSubtraction",
                                       "bias_subtraction_starting",
                                        "get_min_max"),
    
    "get_min_max":                    ("GetMinMax",
                                       "getting_vmin_vmax",
                                       "remove_overscan_regions"),
    
    "remove_overscan_regions":        ("RemoveOverscan",
                                       "removing_overscan_regions",
                                       "create_bottom_row"),
    
    # Creates mosaic
    "create_bottom_row":              ("CreateBottomRow",
                                       "creating_bottom_row",
                                       "create_top_row"),
    
    "create_top_row":                 ("CreateTopRow",
                                       "creating_top_row",
                                       "stack_rows"),
    
    "stack_rows":                     ("StackRows",
                                       "stacking_rows",
                                       "create_mosaic"),
    
    "create_mosaic":                  ("DetectorMosaic",
                                       "creating_full_detector_mosaic",
                                       None),
}