"""
|=============================================================================|
|                                 I N P U T S                                 |
|=============================================================================|
|                                                                             |
| This module contains a routine that reads in the input file and processes   |
| its arguments.                                                              |
|                                                                             |
| Contains                                                                    |
| --------                                                                    |
|     initialise_default_params                                               |
|     write_defaults_to_file                                                  |
|     parse_input                                                             |
|     report_input_file_error                                                 |
|     convert_values                                                          |
|     handle_dependencies                                                     |
|     check_positive_float                                                    |
|     check_list_positive_float                                               |
|     check_positive_int                                                      |
|     check_list_positive_int                                                 |
|     check_int_on_list                                                       |
|     check_unsigned_n_bit_int                                                |
|     check_str_on_list                                                       |
|     check_boolean                                                           |
|                                                                             |
|-----------------------------------------------------------------------------|
| Paul Sharp 26/02/2018, modified by Michail Theofilatos 17/11/2019           |
|=============================================================================|
"""

import ase
import os
import yaml

ase_chemical_symbols = ase.data.chemical_symbols


# =============================================================================
# =============================================================================
def initialise_default_params(calc_name):
    """
    This routine sets up the params dictionary with all of its keywords and
    their default values.

    Returns
    -------
    params : dict
        Dictionary containing each parameter with its default value.

    ---------------------------------------------------------------------------
    Michail Theofilatos 30/10/2019
    """

    # Set up initial parameter dictionary, with all default values.
    params = {

        #
        # General Inputs
        #
        "grid_points": {"value": "1, 1, 1",
                        "specified": False,
                        "description": 'The number of points in each dimension of the grid',
                       },
        "grid_spacing": {"value": "0.5, 0.5, 0.5",
                         "specified": False,
                         "description": 'The spacing between two consecutive grid points on each dimmension.',
                        },
        "angles": {"value": "90,90",
                        "specified": False,
                        "description": 'The minimum and maximum angle between the facets. Default: Orthorhombic',
                        },
        "temp": {"value": 0.0,
                 "specified": False,
                 "description": 'The Monte-Carlo temperature (strictly, the value of kT in eV). Determines whether swaps to basins of higher energy are accepted. Default: 0.0',
                },
        "operations": {"value": 3,
                       "specified": False,
                       "description": 'The operations to be executed in each phase. Default: 3',
                      },
        "operation_1": {"value": "random_structure",
                       "specified": False,
                       "description": 'The first operation to be executed.',
                       },
        "operation_2": {"value": "axes",
                        "specified": False,
                        "description": 'The second operation to be executed.',
                       },
        "operation_3": {"value": "relax_structure",
                        "specified": False,
                        "description": 'The third operation to be executed.',
                       },
        "operation_4": {"value": "relax_structure",
                        "specified": False,
                        "description": 'The fourth operation to be executed.',
                        },
        "operation_5": {"value": "relax_structure",
                        "specified": False,
                        "description": 'The fifth operation to be executed.',
                        },
        "calculator": {"value": "gulp",
                       "specified": False,
                       "description": 'The materials modelling code used for calculations. Default: gulp',
                      },
        "phases": {"value": 1,
                   "specified": False,
                   "description": 'The number of phases to be executed. In each phase, all operations, in the specified order, are executed.',
                  },
        "initial_structure_file": {"value": "initial_structure.cif",
                   "specified": False,
                   "description": 'If specified, read in the initial structure from this cif file.',
                  },
        "initial_structure_folder": {"value": "folder/to/initial/structures",
                   "specified": False,
                   "description": 'If specified, read in the initial structures from this folder.',
                   },
        "function": {"value": "optimization",
                   "specified": False,
                   "description": 'Choice of what to execute.',
                    },
        "global_minimum_energy": {"value": -10.0,
                     "specified": False,
                     "description": 'This parameter specifies the global minimum energy and is used in the "compare_basin_hopping_axes" function.',
                     },
        "whiskers_coarse": {"value": "-15, -20, -25, -30",
                         "specified": False,
                         "description": 'The ranges of the "corse" whisker plot in the "compare_basin_hopping_axes" and "read_results_from_file" functions.',
                         },
        "whiskers_fine": {"value": "-30.0, -30.8, -31.4, -31.6, -31.69",
                            "specified": False,
                            "description": 'The ranges of the "fine" whisker plot in the "compare_basin_hopping_axes" and "read_results_from_file" functions.',
                            },
        "pkl_file": {"value": "C://",
                          "specified": False,
                          "description": 'Path to the .pkl file to plot the data (used in "read_results_from_file" function).',
                          },

        #
        # GULP Inputs
        #
        "calculator_cores": {"value": 1,
                             "specified": False,
                             "description": 'The number of parallel cores used for the calculator. Default: 1.',
                            },
        "gulp_executable": {"value": "/gulp/Src/gulp",
                            "specified": False,
                            "description": 'The filepath of the GULP executable to be used. Default: "./gulp/Src/gulp".',
                           },
        "gulp_library": {"value": "",
                         "specified": False,
                         "description": 'Library file containing the forcefield to be used in GULP calculations. NOTE -- this takes precedence over a library specified in "gulp_options".',
                        },
        "calculator_time_limit": {"value": 10,
                                  "specified": False,
                                  "description": 'Used in the bash "timeout" command. GULP calculations will automatically terminate after this amount of time has expired.',
                                 },
        "gulp_keywords_calculate_energy": {"value": "single, buckingham",
                                           "specified": False,
                                           "description": 'Comma-separated list of keywords for calculating the energy using GULP. Default: "single, buckingham"',
                                          },
        "gulp_keywords_relax": {"value": "opti, pot",
                                "specified": False,
                                "description": 'Comma-separated list of keywords for all GULP relaxations. Default: "opti, pot"',
                                },
        "gulp_options_relax": {"value": [],
                               "specified": False,
                               "description": 'Options for relaxing structure using GULP. Default: None',
                              },
        "gulp_keywords_relax_unit_cell": {"value": "opti, pot, cellonly",
                                          "specified": False,
                                          "description": 'Comma-separated list of keywords for relaxing the unit cell only using GULP. Default: "opti, pot, cellonly"',
                                          },
        "gulp_options_relax_unit_cell": {"value": [],
                                         "specified": False,
                                         "description": 'Options for relaxing unit cell using GULP. Default: None',
                                         },
        "gulp_max_gnorm": {"value": "",
                           "specified": False,
                           "type": "float specified",
                           "description": 'If specified, terminate a GULP calculation if the final gnorm exceeds this value after the first stage.',
                          },

        #
        # Other Inputs
        #
        "atoms_file": {"value": calc_name + ".atoms",
                       "specified": False,
                       "description": 'File in which the species, number, oxidation state and radius of the atoms used in this calculation are specified.',
                      },
        "output_file": {"value": "Results/" + calc_name.rsplit('/', 1)[len(calc_name.rsplit('/', 1))-1] + ".cnplus",
                        "specified": False,
                        "description": 'Output file for this calculation.',
                       },
        "output_trajectory": {"value": "False",
                              "specified": False,
                              "description": 'If true, write ASE trajectory files. Default: True',
                             },
        "seed_bits": {"value": 64,
                      "specified": False,
                      "description": 'The number of bits used in the seed of the random number generator. The allowed values are 32 and 64. Default: 64',
                     },
        "random_seed": {"value": 42,
                        "specified": False,
                        "description": 'The value used to seed the random number generator. Alternatively, the code can generate one itself, which is the default behaviour.',
                       },
        "rng_warm_up": {"value": 10,
                        "specified": False,
                        "description": 'Number of values from the RNG to generate and discard after seeding the generator. Default: 0.',
                       },

        "cp_stacking_sequence": {"value": "",
                                 "specified": False,
                                 "description": 'Anion layer stacking sequence for close packed grids.',
                                },
        "cp_cell_setting": {"value": "primitive",
                            "specified": False,
                            "description": 'Cell setting for close packed grids. Supported values are: "primitive" (default) and "orthorhombic"',
                           },

        "best_struct_file": {"value": "Results/best_structure.gin",
                             "specified": False,
                             "description": 'GULP File to which we write the lowest-energy structure found throughout the run. Default: "best_structure.gin".',
                            }
    }

    return params


# =============================================================================
def write_defaults_to_file(input_file):
    """
    This routine writes the default values for all code options to an input file.

    Parameters
    ----------
    input_file : string
        Name of input file we will write the default values to.

    Returns
    -------
    None

    ---------------------------------------------------------------------------
    Paul Sharp 24/05/2017
    """

    calc_name = input_file.split(".")[0]
    params = initialise_default_params(calc_name)

    # Write to file
    with open(input_file, mode="w") as new_input_file:
        for key in sorted(params.keys()):
            new_input_file.write('%s=%s\n' % (key, params[key]["value"]))

    return None


# =============================================================================
def parse_input(input_file, calc_name):
    """
    This routine runs through the input file for this run of the code, checks
    all the keywords are valid, and applies all values.

    Parameters
    ----------
    input_file : string
        Name of input file for this run of the code.
    calc_name : string
        Name of this job -- default value for name of atoms file.

    Returns
    -------
    params : dict
        Dictionary containing each parameter with its value.
    errors : string
        List of error messages for invalid options.

    ---------------------------------------------------------------------------
    Paul Sharp 28/09/2017
    """

    # Set up initial parameter dictionary, with all default values.
    params = initialise_default_params(calc_name)
    errors = []

    # Read in input file
    with open(input_file, mode="r") as in_file:
        file_contents = in_file.readlines()

    # Remove entries that start with a comment character -- either # or ! -- or newline
    file_contents[:] = [x for x in file_contents if not x.startswith(("#", "!", "\n"))]

    # Remove spaces from each entry and newlines -- exclude GULP options because they will contain spaces
    for i in range(0, len(file_contents)):
        if ("gulp" and "options") not in file_contents[i]:
            file_contents[i] = ''.join(file_contents[i].split())

    # Check whether the keywords are valid
    invalid_keywords = []
    for entry in file_contents[:]:
        keyword = entry.split("=")[0].strip()

        if keyword in params:
            params[keyword]["specified"] = True
        else:
            invalid_keywords.append(keyword)
            file_contents.remove(entry)

    if len(invalid_keywords) > 0:
        errors.append("The keyword(s) %s are invalid."
                      % (", ".join(invalid_keywords)))

    # Apply the new values to the keywords
    for entry in file_contents:

        keyword = entry.split("=")[0].strip()
        value = entry.split("=")[1].strip()

        params[keyword]["value"] = value

    return params, errors


# =============================================================================
def report_input_file_errors(errors, input_file, error_file):
    """
    This routine prints the error list for an input file in the correct format.

    Parameters
    ----------
    errors : string
        List of error messages for invalid options.
    input_file : string
        Name of input file that was parsed.
    error_file : string
        Name of the file in which we report the errors.

    Returns
    -------
    None

    ---------------------------------------------------------------------------
    Paul Sharp 07/09/2017
    """

    with open(error_file, mode="a") as error_log:

        error_log.write('%i errors were found in the input file: "%s"\n'
                        % (len(errors), input_file))

        for error in errors:
            error_log.write('\tERROR in input file "%s" -- %s\n'
                            % (input_file, error))

        error_log.write('\n')

    return None


# =============================================================================
def convert_values(params, errors):
    """
    This routine runs through the dictionary of input parameters, converts them
    to the correct type, and ensures they are valid.

    Parameters
    ----------
    params : dict
        Dictionary containing each parameter with its value.
    errors : string
        List of error messages for invalid options.

    Returns
    -------
    params : dict
        Dictionary containing each parameter with its value.
    errors : string
        Updated list of error messages for invalid options.

    ---------------------------------------------------------------------------
    Paul Sharp 26/02/2018
    """

    supported_calculators = ["gulp"]

    # Lists
    params["grid_points"]["value"] = params["grid_points"]["value"].split(",")
    params["angles"]["value"] = params["angles"]["value"].split(",")
    params["grid_spacing"]["value"] = params["grid_spacing"]["value"].split(",")
    params["whiskers_coarse"]["value"] = params["whiskers_coarse"]["value"].split(",")
    params["whiskers_fine"]["value"] = params["whiskers_fine"]["value"].split(",")


    # Booleans
    boolean_keys = ["output_trajectory"]

    for key in boolean_keys:
        params[key]["value"], errors = check_boolean(key, params[key]["value"], errors)

    # Integers
    integer_keys = ["operations", "phases", "calculator_cores", "calculator_time_limit", "grid_points", "angles", "rng_warm_up"]

    for key in integer_keys:
        params[key]["value"], errors = check_positive_int(key, params[key]["value"], errors)

    # Floats
    positive_float_keys = ["grid_spacing", "temp"]
    other_float_keys = ["whiskers_coarse", "whiskers_fine", "global_minimum_energy"]

    for key in positive_float_keys:
        params[key]["value"], errors = check_positive_float(key, params[key]["value"], errors)
    for key in other_float_keys:
        params[key]["value"], errors = check_float(key, params[key]["value"], errors)


    # Strings
    # Easy test for existence of executable would be nice -- only exists in python 3
    string_keys = ["operation_1", "operation_2", "operation_3", "operation_4", "operation_5", "initial_structure_file",
                   "initial_structure_folder", "gulp_keywords_calculate_energy", "gulp_keywords_relax", "gulp_keywords_relax_unit_cell",
                   "atoms_file", "best_struct_file", "cp_cell_setting",
                   "gulp_executable", "gulp_library", "output_file", "function"]

    for key in string_keys:
        params[key]["value"] = str(params[key]["value"])

    params["calculator"]["value"], errors, calculator_valid = check_str_on_list("calculator", params["calculator"]["value"], errors, supported_calculators)
    #params["cp_cell_setting"]["value"], errors, cell_setting_valid = check_str_on_list("cp_cell_setting", params["cp_cell_setting"]["value"], errors, supported_cell_settings)


    # Lists
    # Let GULP deal with full verification of it's keywords and options
    list_to_string_keys = ["gulp_keywords_relax", "gulp_keywords_relax_unit_cell", "gulp_keywords_calculate_energy"]

    for key in list_to_string_keys:
        params[key]["value"] = params[key]["value"].replace(",", " ")

    list_keys = ["gulp_options_relax", "gulp_options_relax_unit_cell"]

    for key in list_keys:
        try:
            params[key]["value"] = [option.strip() for option in params[key]["value"].split(",")]
        except AttributeError:
            pass

    return params, errors


# =============================================================================
def handle_dependencies(params, errors):
    """
    This routine runs through the dictionary of input parameters,
    cross-checking options that depend on others.

    Parameters
    ----------
    params : dict
        Dictionary containing each parameter with its value.
    errors : string
        List of error messages for invalid options.

    Returns
    -------
    params : dict
        Dictionary containing each parameter with its value.
    errors : string
        Updated list of error messages for invalid options.

    ---------------------------------------------------------------------------
    Paul Sharp 06/11/2017
    """

    supported_seed_bits = [32, 64]

    # Random seed must be the correct number of bits.
    params["seed_bits"]["value"], errors, seed_bits_valid = check_int_on_list("seed_bits", params["seed_bits"]["value"], errors, supported_seed_bits)
    if seed_bits_valid:
        params["random_seed"]["value"], errors = check_unsigned_n_bit_int("random_seed", params["random_seed"]["value"], errors, params["seed_bits"]["value"])
    else:
        errors.append('"random_seed" has not been verified because it must be either a 32- or 64-bit integer, but an invalid number of bits was entered.')

    # We need an ".atoms" file -- check the file given is valid and exists
    if (params["atoms_file"]["value"].split(".")[1] != "atoms"):
        errors.append('the ".atoms" file is specified as "%s", but that file does not have the ".atoms" extension.'
                      % (params["atoms_file"]["value"]))
    if (not os.path.isfile(params["atoms_file"]["value"])):
        errors.append('the ".atoms" file is specified as "%s", but that file does not exist.'
                      % (params["atoms_file"]["value"]))

    # If we are initialising from a ".cif" file -- check the file given is valid and exists
    if params["initial_structure_file"]["specified"]:
        if (params["initial_structure_file"]["value"].split(".")[1] != "cif"):
            errors.append('the ".cif" file for the initial structure is specified as "%s", but that file does not have the ".cif" extension.'
                          % (params["initial_structure_file"]["value"]))
        if (not os.path.isfile(params["initial_structure_file"]["value"])):
            errors.append('the ".cif" file for the initial structure is specified as "%s", but that file does not exist.'
                          % (params["initial_structure_file"]["value"]))

    if params["initial_structure_folder"]["specified"]:
        if (not os.path.isdir(params["initial_structure_folder"]["value"])):
            errors.append('the initial structures folder ("%s") does not exist.'
                          % (params["initial_structure_folder"]["value"]))

    return params, errors


# =============================================================================
def check_positive_float(keyword, value, errors):
    try:
        value = float(value)
    except ValueError:
        errors.append('"%s" is %s, but should be a positive number.'
                      % (keyword, value))
    # A TypeError means that a list was inputted
    except TypeError:
        value, errors = check_list_positive_float(keyword, value, errors)
    else:
        if value < 0.0:
            errors.append('"%s" is %s, but should be a positive number.'
                          % (keyword, value))

    return value, errors


# =============================================================================
def check_list_positive_float(keyword, value, errors):
    try:
        value = [float(x) for x in value]
    except ValueError:
        errors.append('"%s" is %s, but all values should be positive floats.'
                      % (keyword, value))
    else:
        if any(i < 0.0 for i in value):
            errors.append('"%s" is %s, but all values should be positive floats.'
                          % (keyword, value))

    return value, errors

# =============================================================================

def check_float(keyword, value, errors):
    try:
        value = float(value)
    except ValueError:
        errors.append('"%s" is %s, but should be a positive number.'
                      % (keyword, value))
    # A TypeError means that a list was inputted
    except TypeError:
        value, errors = check_list_float(keyword, value, errors)

    return value, errors

# =============================================================================

def check_list_float(keyword, value, errors):
    try:
        value = [float(x) for x in value]
    except ValueError:
        errors.append('"%s" is %s, but all values should be positive floats.'
                      % (keyword, value))
    return value, errors


# =============================================================================
def check_positive_int(keyword, value, errors):
    try:
        value = int(value)
    except ValueError:
        errors.append('"%s" is %s, but should be a positive integer.'
                      % (keyword, value))
    # A TypeError means that a list was inputted
    except TypeError:
        value, errors = check_list_positive_int(keyword, value, errors)
    else:
        if value < 0:
            errors.append('"%s" is %s, but should be a positive integer.'
                          % (keyword, value))

    return value, errors


# =============================================================================
def check_list_positive_int(keyword, value, errors):
    try:
        value = [int(x) for x in value]

    except ValueError:
        errors.append('"%s" is %s, but all values should be positive integers.'
                      % (keyword, value))
    else:
        if any(i < 0 for i in value):
            errors.append('"%s" is %s, but all values should be positive integers.'
                          % (keyword, value))

    return value, errors


# =============================================================================
def check_int_on_list(keyword, value, errors, allowed_ints):
    valid = True
    try:
        value = int(value)
    except ValueError:
        errors.append('"%s" is %s, but should be one of the supported values: %s'
                      % (keyword, value, ', '.join([str(x) for x in allowed_ints])))
        valid = False
    else:
        if value not in allowed_ints:
            errors.append('"%s" is %s, but should be one of the supported values: %s'
                          % (keyword, value, ', '.join([str(x) for x in allowed_ints])))
            valid = False

    return value, errors, valid


# =============================================================================
def check_unsigned_n_bit_int(keyword, value, errors, bits):
    try:
        value = int(value)
    except ValueError:
        errors.append('"%s" is %s, but should be an unsigned %i-bit integer.'
                      % (keyword, value, bits))
    else:
        if (value < 0) or (value > (2 ** int(bits)) - 1):
            errors.append('"%s" is %s, but should be an unsigned %i-bit integer.'
                          % (keyword, value, bits))

    return value, errors


# =============================================================================
def check_str_on_list(keyword, value, errors, allowed_strings):
    valid = True
    if value not in allowed_strings:
        errors.append('"%s" is %s, but should only be one of the supported values: %s'
                      % (keyword, value, ', '.join([str(x) for x in allowed_strings])))
        valid = False

    return value, errors, valid


# =============================================================================
def check_boolean(keyword, value, errors):
    true_synonyms = ["true", "t", "yes", "on"]
    false_synonyms = ["false", "f", "no", "off"]

    if value.lower() in true_synonyms:
        value = True
    elif value.lower() in false_synonyms:
        value = False
    else:
        errors.append('"%s" is %s, but should be "True" or "False".'
                      % (keyword, value))

    return value, errors

