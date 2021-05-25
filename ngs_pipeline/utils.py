def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values
    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False



def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "None", "none", "F", "f"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False



def get_chromsizes_from_ucsc(
    genome,
    saveas=None,
    mysql="mysql",
    fetchchromsizes="fetchChromSizes",
    timeout=None,
    host_url="genome-mysql.cse.ucsc.edu",
):
    """
    Download chrom size info for *genome* from UCSC and returns the dictionary.

    Taken from pybedtools.

    Parameters
    ----------

    genome : str
        Name of the genome assembly (e.g., "hg38")

    saveas : str
        Filename to save output to. Dictionary will still be returned.

    mysql, fetchchromsizes : str
        Paths to MySQL and fetchChromSizes.

    timeout : float
        How long to wait for a response; mostly used for testing.

    host_url : str
        URL of UCSC mirror MySQL server.
    """
    if not internet_on(timeout=timeout):
        raise ValueError(
            "It appears you don't have an internet connection "
            "-- unable to get chromsizes from UCSC"
        )
    cmds = [
        mysql,
        "--user=genome",
        "--host=" + host_url,
        "-A",
        "-e",
        "select chrom, size from %s.chromInfo" % genome,
    ]
    failures = []
    d = {}
    try:
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=BUFSIZE
        )
        stdout, stderr = p.communicate()
        if stderr:
            print(stderr)
            print("Commands were:\n")
            print((subprocess.list2cmdline(cmds)))

        lines = stdout.splitlines()[1:]
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode("UTF-8")
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

    except OSError as err:
        if err.errno == 2:
            failures.append("Can't find mysql at path {0}".format(mysql))
        else:
            raise
    try:
        cmds = [fetchchromsizes, genome]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=BUFSIZE
        )
        stdout, stderr = p.communicate()
        if stderr:
            if "INFO: trying WGET" not in str(stderr):
                print(stderr)
                print("Commands were:\n")
                print((subprocess.list2cmdline(cmds)))

        lines = stdout.splitlines()
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode("UTF-8")
            chrom, size = line.split()
            d[chrom] = (0, int(size))

        if saveas is not None:
            chromsizes_to_file(d, saveas)

    except OSError as err:
        if err.errno == 2:
            failures.append("Can't find path to fetchChromsizes")

    if not d:
        raise OSError(failures)
    return d


