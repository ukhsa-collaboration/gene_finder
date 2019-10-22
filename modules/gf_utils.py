import logging
import sys
import os
import time

# --------------------------------------------------------------------------- #

class PheException(Exception):
    '''
    This is the top level class that EVERYTHING must be derived from. In particular,
    this class contains an abstract property called 'phe_return_code'. This property
    must be implemented and the individual implementation will have it's own
    exit code. which will be propogated to the calling functions, if needs be.

    PheException must not be passed as is.
    '''

    def __init__(self, msg, cause, phe_return_code=255):
        '''
        Constructor
        '''
        super(Exception, self).__init__(msg)
        self._phe_return_code = phe_return_code
        self._cause = cause

    @property
    def phe_return_code(self):
        '''
        Read-only attribute that holds the return status that should be exited with.
        '''
        return self._phe_return_code

    @property
    def cause(self):
        '''
        Read-only attribute that indicates the root cause of the exception raised.
        '''
        return self._cause

# --------------------------------------------------------------------------- #

class PheExternalError(PheException):
    '''
    Exception class designed to be raised when an external command/process
        fails. Instead of falling over quietly, this exception can be raised. The
        exception includes the message to be put into the logs and the cause of
        the exception. In this case, the cause should generally be subprocess.CallerProcessError.
        The particulars of the failed command can be found inside the cause.
        If the catcher of this exception choses to exit the code, 'phe_return_code'
        should be used to indicate the cause of it all.
    '''
    def __init__(self, msg, cause):
        '''
        Constructor for the PheExternalError

        @param msg: Message to be displayed with the exception.
        @type msg: str.
        @param cause: Cause of this exception, usually subprocess.CalledProcessError.
        @type cause: class.
        '''
        super(PheExternalError, self).__init__(msg, cause, 55)

# --------------------------------------------------------------------------- #

def getNGSSproductionHandlers(logs_dir, component_name):
    '''
    returns logging handlers conferring on loggers to which they are added
    the logging behaviours expected within Production NGSService pipelines
    args
    ----
    logs_dir:str, full path with 'logs' endpoint
    component_name: str
    returns
    -------
    tuple of handlers writing
        log.debug to console
        log.info, log.warning to stdout log
        log.error, log.exception, log.critical to stderr log
    # note that the log.warning write into stdout log differs
       from the setup_logger() function which directs warnings
       into the stderr log
    '''
    class LogLevelFilter(object):
        def __init__(self, level):
            self.__level = level
        def filter(self, logRecord):
            return logRecord.levelno <= self.__level

    ngs_formatter = logging.Formatter('%(asctime)s\n%(message)s')

    error_file = os.path.sep.join([logs_dir, component_name+'.stderr'])
    fh_stderr = logging.FileHandler(error_file)
    fh_stderr.setFormatter(ngs_formatter)
    fh_stderr.setLevel(logging.ERROR)
    fh_stderr.addFilter(LogLevelFilter(logging.ERROR))

    info_file = os.path.sep.join([logs_dir, component_name+'.stdout'])
    fh_stdout = logging.FileHandler(info_file)
    fh_stdout.setFormatter(ngs_formatter)
    fh_stdout.setLevel(logging.INFO)
    fh_stdout.addFilter(LogLevelFilter(logging.INFO))

    lvl_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    fh_warn = logging.FileHandler(info_file)
    fh_warn.setFormatter(lvl_formatter)
    fh_warn.setLevel(logging.WARNING)
    fh_warn.addFilter(LogLevelFilter(logging.WARNING))

    fh_crit = logging.FileHandler(error_file)
    fh_crit.setFormatter(lvl_formatter)
    fh_crit.setLevel(logging.CRITICAL)
    fh_crit.addFilter(LogLevelFilter(logging.CRITICAL))

    sh_formatter = \
    logging.Formatter('::%(levelname)s::%(asctime)s %(name)s %(message)s')
    sh_debug = logging.StreamHandler(sys.stdout)
    sh_debug.setFormatter(sh_formatter)
    sh_debug.setLevel(logging.DEBUG)
    sh_debug.addFilter(LogLevelFilter(logging.DEBUG))

    return fh_stderr, fh_stdout, fh_warn, fh_crit, sh_debug

# --------------------------------------------------------------------------- #

def setup_logging(logs_dirname, time_stamp=None, project_basename=None, \
                    logger_name=None, loglevel=None):
    '''
    Wrapper to setup_logger which abstracts the details &
    provides for time-stamping.

    Parameters
    ----------
    logs_dirname : str
        Path to which logs/ dir will be appended
    time_stamp : None or str
        If None, a stamp-stamp in the format 'yyyy-mm-dd-hm' is generated,
        otherwise the timestamp passed is used to name the output log files.
    project_basename : None or str
        Basename of the calling project; used to name the logger unless a
        logger_name is passed
    logger_name : None or str
        Default is None, in which case the logger_name is
        project_basename_logger
        Otherwise the passed logger_name is used.

    Return
    ------
    logger : An instance of the logger class

    Note
    ----
    When project_basename is passed the logger_name is project_basename_logger &
    can be recovered elsewhere in the project using the statement
        logger = logging.getlogger( logger_name )
    otherwise the passed logger_name is used.

    When project_basename is passed the log files are written to
    log_dirname/logs/ as project_basename_timestamp.stdout &
    project_basename_timestamp.stderr
    otherwise the passed logger_name is used to write to files
    logger_name_timestamp.stdout & logger_name_timestamp.stdout

    '''

    logs_path = os.path.join(logs_dirname, 'logs/')
    #print 'logs_path:', logs_path
    if not os.path.exists(logs_path):
        os.makedirs(logs_path)

    if not time_stamp:
        time_stamp = time.strftime('%Y-%m-%d_%H%M%S')

    if not logger_name:
        logger_name = '_'.join([project_basename, 'logger'])
        stdout, stderr = [''.join([logs_path, project_basename, '_', time_stamp, log]) for log in ['.stdout', '.stderr']]
    else:
        stdout, stderr = [''.join([logs_path, logger_name, '_', time_stamp, log]) for log in ['.stdout', '.stderr']]

    logger = setup_logger(info_file=stdout, error_file=stderr, logger_name=logger_name, loglevel=loglevel)
    return logger


# --------------------------------------------------------------------------- #


def setup_logger(info_file="stdout.log", error_file="stderr.log", \
                logger_name='stdout_stderr_logger', loglevel=None):
    """
    A function to set up a logger for writing out log files

    Parameters
    ----------
    info_file : String
        Path to info level log file (default: stdout.log)
    error_file : String
        Path to error level log file (default: stderr.log)
    logger_name : String
        Name for logger (default: 'stdout_stderr_logger')

    Returns
    -------
    logger : A instance of the Logger class

    Note
    ----
    In a code context into which 'logging' has been imported, the logger_name
    parameter allows a logger instantiated elsewhere to be recovered using
    logger = logging.getLogger( 'stdout_stderr_logger' ), obviating the need
    to pass logger as an arg to any function.

    """

    class LogLevelFilter(object):
        def __init__(self, level):
            self.__level = level

        def filter(self, logRecord):
            return logRecord.levelno <= self.__level

    logger = logging.getLogger(logger_name)
    if not loglevel:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(loglevel)

    formatter = logging.Formatter('%(asctime)s\n%(message)s')

    handler_stderr = logging.FileHandler(error_file)
    handler_stderr.setLevel(logging.ERROR)
    handler_stderr.setFormatter(formatter)
    handler_stderr.addFilter(LogLevelFilter(logging.ERROR))
    logger.addHandler(handler_stderr)

    handler_stdout = logging.FileHandler(info_file)
    handler_stdout.setLevel(logging.INFO)
    handler_stdout.setFormatter(formatter)
    handler_stdout.addFilter(LogLevelFilter(logging.INFO))
    logger.addHandler(handler_stdout)

    # add handlers for other log_levels
    # formatter reporting loglevel before logged msg
    lvl_formatter = logging.Formatter('%(asctime)s\n%(levelname)s::%(message)s')
    for level, output in [(logging.WARNING, error_file),
                          (logging.CRITICAL, error_file)]:
        handler = logging.FileHandler(output)
        handler.setLevel(level)
        handler.setFormatter(lvl_formatter)
        handler.addFilter(LogLevelFilter(level))
        logger.addHandler(handler)
    # stream DEBUG to console
    inline_formatter = logging.Formatter('\n::%(levelname)s::%(message)s')
    for level, output in [(logging.DEBUG, sys.stdout)]:
        handler = logging.StreamHandler(output)
        handler.setLevel(level)
        handler.setFormatter(inline_formatter)
        handler.addFilter(LogLevelFilter(level))
        logger.addHandler(handler)

    return logger

# --------------------------------------------------------------------------- #

def write_log(logger, log_text, log_level):
    """
    Writes text to a logger at a particular log level

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    log_text : String
        The text to be written to the log file
    log_level : String
        The level of logging to which the text should be
            written (either 'info' or 'error')

    """

    if log_level == "error":
        logger.error(log_text)
    elif log_level == "info":
        logger.info(log_text)

# --------------------------------------------------------------------------- #

def log_process(logger, process, log_info_to="info",
                log_error_to="error", limit_logging=0):
    """
    A function to log the output of a subprocess.Popen call

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    process : process pipe
        A process created by subprocess.Popen
    log_info_to: String
        The level at which to log info level logs into (default 'info')
    log_error_to: String
        The level at which to log error level logs into (default 'error')
    limit_logging: Integer
        Limit logging to either stdout (1) or stderr (2), \
                log both if 0 (default 0)

    Returns
    -------
    stdout : String
        The stdout from the process
    stderr : String
        The stderr from the process
    """
    # if the process gives a exit status greater than 1
    #, i.e. an genuine error, then log_error_to 'error'.
    if process.returncode > 0:
        log_error_to = "error"

    stdout = ""
    stderr = ""

    if limit_logging != 2:
        if not process.stdout == None:
            stdout = process.stdout.read()
        if len(stdout) > 0:
            write_log(logger, stdout, log_info_to)
    if limit_logging != 1:
        if not process.stderr == None:
            stderr = process.stderr.read()
        if len(stderr) > 0:
            write_log(logger, stderr, log_error_to)

    return stdout, stderr

# --------------------------------------------------------------------------- #

def write_header_to_log(logger, header_text, log_level):
    """
    A utility function to write some header text bound by asterisks to the log
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    log_level : String
        The level at which to log
    """
    write_log(logger, "******** " + header_text + " ********", log_level)

# --------------------------------------------------------------------------- #

def error_header(logger, header_text):
    """
    utility function; writes header_text bound by asterisks to error log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "error")

# --------------------------------------------------------------------------- #

def info_header(logger, header_text):
    """
    utility function; writes header_text bound by asterisks to the info log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "info")

# --------------------------------------------------------------------------- #

def get_logger_path(args, dir_name='logs'):
    '''
    Returns a path for logger output according to params in the return from
    parser.parse_args(). By default returns output_dir/logs, else returns
    input_dir/logs, else fastqfile_dir/logs

    Args:
        args, Namespace object : the return from parser.parse_args()
        dir_name, string : terminal dir appended to derived path

    Returns:
        output_dir, string : fully-specified output path terminating
        with /dir_name to which logs can be written

    Side effect:
        creates path if necessary

    '''

    log_path = None
    # if an output_dir has been passed check it & use it
    if not args.output_dir is None:
        if not os.path.isdir(args.output_dir):
            print('the output_dir passed (' \
                    + args.output_dir + ') does not exist')
            print('making dir: ' + args.output_dir)
            os.makedirs(args.output_dir)
        log_path = args.output_dir

    # otherwise use any input_dir
    elif not args.input_dir is None:
        if not os.path.isdir(args.input_dir):
            print('ERROR: the input_dir passed (' \
                    + args.input_dir + ') is not valid')
        else:
            log_path = args.input_dir

    # with neither input_dir nor output_dir passed
    # extract a path from fastq files passed
    elif args.fastq_1 and args.fastq_2:
        fqs = [args.fastq_1, args.fastq_2]
        fq_path = set([os.path.split(fq)[0] for fq in fqs])
        if not len(fq_path) == 1:
            print('ERROR: fastq files passed (' \
                    + str(fqs) + ') have different paths')
        else:
            (fq_path,) = fq_path
            if not os.path.isdir(fq_path):
                print(''.join(['fastq files passed have invalid path: ', \
                        fq_path]))
            else:
                log_path = fq_path

    assert log_path is not None, \
            'a valid path for logger could not be derived from the passed args'

    full_log_path = ''.join([log_path, '/', dir_name])
    if not os.path.isdir(full_log_path):
        os.mkdir(full_log_path)
    return full_log_path


# ##############################################################################################

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
        (a useful default = logger.handlers[0].baseFilename, which returns
            the stderr.log FileHandler added first by log_writer.setup_logger)
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
    import traceback

    try:
        return function(*parameters, **named_parameters)
    except PheException as phe_e:
        # This exception is created by the 'call_external', when exit code != 0

        logger = logging.getLogger("stdout_stderr_logger")
        logger.exception(function.__name__ + " has raised an exception: \n" + str(phe_e))

        # Exit with the returncode specified in the called process.
        # TODO: Should it be a different return code? E.g. ranged for traceback.
        sys.exit(phe_e.phe_return_code)
    except Exception:
        error_string = "There was an error in the function '" + function.__name__ + "'"
        error_divider = "_" * 60
        print (error_string)
        print (error_divider)
        traceback.print_exc()
        print (error_divider)

        error_file = open(error_filepath, "a")
        error_file.write(error_string + "\n")
        error_file.write(error_divider + "\n")
        traceback.print_exc(file=error_file)
        error_file.write(error_divider + "\n")
        error_file.close()
        sys.exit(1)
