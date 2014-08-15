def try_and_except(error_filepath, function, *parameters, **named_parameters):
    """
    This wraps a function in try and except clause. If an error is caught this will trigger
    1) reporting of the error to stdout
    2) writing of the error into an error file
    3) a sys.exit with code 1

    Parameters
    ----------
    error_file : String
        path to log file to capture error in
    function: Function
        the function
    parameters : all non-named parameters for the function
    named_parameters : all named paramaters for the function

    Notes
    -----
    Returns the returns from the function called

    Examples
    --------

    assuming a function

    def my_func(a,b,c = None)
        ......
        return d

    This function would be called as follows:

    return_val = try_and_except("stderr.log", my_func, 1, 2, c = 3)
    """
    import sys, traceback
    try:
        return function(*parameters, **named_parameters)
    except Exception as e:
        error_string = "There was an error in the function '" + function.__name__ + "'"
        error_divider = "_" * 60
        print error_string
        print error_divider
        traceback.print_exc()
        print error_divider

        error_file = open(error_filepath, "a")
        error_file.write(error_string + "\n")
        error_file.write(error_divider + "\n")
        traceback.print_exc(file = error_file)
        error_file.write(error_divider  + "\n")
        error_file.close()
        sys.exit(1)

def check_file_exists(filepath, file_description):
    """
    A function to check if a file exists.
    It will print out an error message and exit if the file is not found

    Parameters
    ----------
    filepath : String
        the path to the file to be checked
    file_description : String
        a description of the file to be checked e.g "config file"
    """
    import os
    if not os.path.exists(filepath):
        print("The " + file_description + " (" + filepath + ") does not exist")
        sys.exit(1)
