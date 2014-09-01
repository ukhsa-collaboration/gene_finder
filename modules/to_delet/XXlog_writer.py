import logging

def setup_logger(info_file = "stdout.log", error_file = "stderr.log"):
	"""
	A function to set up a logger for writing out log files

	Parameters
	----------
	info_file : String
		Path to info level log file (default: stdout.log)
	error_file : String
		Path to error level log file (default: stderr.log)

	Returns
	-------
	logger : A instance of the Logger class
	"""

	class LogLevelFilter(object):
	    def __init__(self, level):
	        self.__level = level

	    def filter(self, logRecord):
	        return logRecord.levelno <= self.__level

	logger = logging.getLogger('error_testing')
	logger.setLevel(logging.DEBUG)

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

	return logger

def write_log(logger, log_text,log_level):
	"""
	Writes text to a logger at a particular log level

	Parameters
	----------
	logger : Logger
		An instance of the Logger class
	log_text : String
		The text to be written to the log file
	log_level : String
		The level of logging to which the text should be written (either 'info' or 'error')

	"""

	if log_level == "error":
		logger.error(log_text)
	elif log_level == "info":
		logger.info(log_text)

def log_process(logger, process,log_info_to = "info", log_error_to = "error"):
	"""
	A function to log the output of a subprocess.Popen call

	Parameters
	----------
	logger : Logger
		An instance of the Logger class
	process : process pipe
		A process created by subprocess.Popen
	log_info_to: String
		The level at which to log info level logs ino (default 'info')
	log_error_to: String
		The level at which to log error level logs ino (default 'error')

	Returns
	-------
	stdout : String
		The stdout from the process
	stderr : String
		The stderr from the process
	"""
	#if the process gives a exit status greater than 1, i.e. an genuine error, then log_error_to 'error'.
	if process.returncode > 0:
		log_error_to = "error"
	if not process.stdout == None:
		stdout = process.stdout.read()
	else:
		stdout = ""
	if not process.stdout == None:
		stderr = process.stderr.read()
	else:
		stderr = ""
	if len(stdout) > 0:
		write_log(logger, stdout, log_info_to)
	if len(stderr) > 0:
		write_log(logger, stderr, log_error_to)
	return stdout, stderr

def info_header(logger, header_text):
	"""
	A utility function to write some header text bound by asterisks to the info log_text
	Parameters
	----------
	logger : Logger
		An instance of the Logger class
	header_text : String
		The header text to embed within asterisks
	"""
	write_log(logger, "******** " + header_text + " ********", "info")