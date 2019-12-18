
import collections
import re

def read_table (lines):
    """Reads to the next divider line.
    
    Read lines are returned, 'lines' is reduced to everything remaining.
    """
    table_lines = []
    while len(lines):
        popline = lines[0]
        lines = lines[1:]
        if popline[:6] == "======":
            break
        table_lines += [popline]
    return table_lines


def identify_tables (lines):
    """Splits the input into a list of tables."
    
    Tables are identified by a divider line followed by a line like
    
        Table xxx : The Table Title
    
    They end at the next divider line.
    """
    tables = {}
    while len(lines):
        popline = lines[0]
        lines = lines[1:]
        if popline[:6] == "======" and len(lines):
            popline = lines[0]
            lines = lines[1:]
            m = re.search ('\ATable\s.*\:\s*(\w.*\w)\Z', popline)
            if m:
                table_name = m.group(1)
                tables[table_name] = read_table(lines)
    return tables


def subtables (table):
    subtables_list = []
    # Skip to first subtable
    i = 0
    while i < len(table):
        if table[i][:6] == "------":
            i += 1
            break
        i += 1
    while i < len(table):
        subtables_list += [[]]
        while i < len(table):
            if table[i][:6] == "------":
                i += 1
                break
            subtables_list[-1] += [table[i]]
            i += 1
    return subtables_list


def subtable_by_key (table, key, value):
    for st in subtables (table):
        for l in st:
            columns = l.split(':')
            if len(columns) == 2:
                if columns[0].strip() == key and columns[1].strip() == value:
                    return st
    return None
        


def get_labelled_value (label, table):
    """Looks for lines with the format
    
        label:   value
        
    and returns value.
    """
    for l in table:
        m = re.search ('%s\s*:?\s*((?:[\-\w].*\w)|([\.\\\\/\w]+))' % label, l)
        if m:
            return m.group(1)
    return None


def get_labelled_triplet (label, table):
    """Looks for lines with the format
    
        label:   value value value
        
    and returns value, value, value.
    """
    for l in table:
        m = re.search ('%s\s*:?\s*([\-\+\.01-9E]+|NAN|INF)\s+([\-\+\.01-9E]+|NAN|INF)\s+([\-\+\.01-9E]+|NAN|INF)' % label, l)
        if m:
            return (m.group(1), m.group(2), m.group(3))
    return None


def read_simple_table_by_row (first_header, row, column, table):
    """Reads a simple table that starts with a header.
    
    Once the header is found, the value in the specified row and column
    will be returned.
    """
    # First find header
    found_header = False
    for i in xrange(len(table)):
        if table[i].split()[0].lower() == first_header.lower():
            found_header = True
            break
    if not found_header:
        return None
    return table[i+row].split()[column]


def read_simple_table_by_key (first_header, row_key, column, table):
    """Reads a simple table that starts with a header.

    Once the header is found, the value in the specified row and column
    will be returned.
    """
    # First find header
    found_header = False
    for i in xrange(len(table)):
        try:
            if table[i].split()[0].lower() == first_header.lower():
                found_header = True
                break
        except:
            pass
    if not found_header:
        return None
    for i in xrange(i,len(table)):
        columns = table[i].split()
        if columns[0] == row_key:
            return columns[column]
    return None


def alternates (functions):
    for f in functions:
        try:
            value = f()
            if value != None:
                return value
        except:
            pass
    raise Exception("No alternate look-up resolved for variable.")

tables = {}
def get_tables (input_file):
    """Opens the specified file, strips lines of white space, and
        calls identify_tables."""
    global tables
    
    lines = []
    for l in open(input_file):
        # Strip on right and at most one space on left
        if len(l) > 0 and l[0] == " ":
            lines += [l[1:].rstrip()]
        else:
            lines += [l.rstrip()]
    tables = identify_tables (lines)
    return tables

"""Define a dictionary of functions to use to find values.
"""
lookups = collections.OrderedDict()
lookups['filename'] = \
        lambda : get_labelled_value ("Filename", tables["Model Input"])
lookups["dx_els"] = \
    lambda : get_labelled_value ("Element dim X", tables["Model Input"])
lookups["dy_els"] = \
    lambda : get_labelled_value ("Element dim Y", tables["Model Input"])
lookups["dz_els"] = \
    lambda : get_labelled_value ("Element dim Z", tables["Model Input"])
lookups["num_els"] = \
        lambda : get_labelled_value ("Number of elements", tables["Model Input"])
lookups["num_nodes"] = \
        lambda : get_labelled_value ("Number of nodes", tables["Model Input"])
lookups["num_nodes_per_els"] = \
        lambda : get_labelled_value ("Number of nodes per element", tables["Model Input"])
lookups["dim_of_problem"] = \
        lambda : get_labelled_value ("Dimension of problem", tables["Model Input"])
lookups["num_mats"] = \
        lambda : alternates (
            [ lambda: get_labelled_value ("Number of materials",
                                            tables["Materials"]),
                lambda: get_labelled_value ("Number of materials available",
                                            tables["Material Table"]) ])
lookups["id_mat%m"] = \
        lambda m : alternates (
            [ lambda : read_simple_table_by_row ("m", m, 1,
                                                tables["Materials"]),
            lambda : read_simple_table_by_row ("id", m, 0,
                                                tables["Material Table"]) ])
lookups["count_mat%m"] = \
        lambda m : alternates (
            [ lambda : read_simple_table_by_row ("m", m, 5,
                                                tables["Materials"]),
            lambda : read_simple_table_by_row ("material", m, 1,
                                                tables["Material Distribution"]) ])
lookups["num_pp_sets"] = \
        lambda : alternates ( 
            [ lambda : get_labelled_value ("Number of sets",
                                            tables["Post-processing Sets"]),
            lambda : get_labelled_value ("Number of node sets",
                                            tables["Nodal Displacements"]) ])
lookups["eps_ave_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_stddev_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_min_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_max_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_skew_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_kurt_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_median_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc05_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc25_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc75_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc95_xx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 1,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_ave_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_stddev_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_min_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_max_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_skew_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_kurt_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_median_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc05_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc25_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc75_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc95_yy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 2,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_ave_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_stddev_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_min_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_max_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_skew_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_kurt_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_median_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc05_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc25_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc75_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["eps_perc95_zz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 3,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_ave_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_stddev_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_min_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_max_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_skew_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_kurt_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_median_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc05_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc25_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc75_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc95_yz"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 4,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_ave_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_stddev_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_min_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_max_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_skew_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_kurt_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_median_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc05_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc25_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc75_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc95_zx"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 5,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_ave_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "average", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_stddev_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "std_dev", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_min_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "minimum", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_max_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "maximum", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_skew_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "skewness", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_kurt_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "kurtosis", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_median_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "median", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc05_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc05", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc25_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc25", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc75_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc75", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["gam_perc95_xy"] = \
        lambda m : read_simple_table_by_key ("epsilon_xx", "perc95", 6,
                        subtable_by_key(tables["Strain"],"m",str(m)))
lookups["sig_ave_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_xx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 1,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_ave_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_yy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 2,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_ave_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_zz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 3,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_ave_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_yz"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 4,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_ave_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_zx"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 5,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_ave_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "average", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_stddev_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "std_dev", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_min_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "minimum", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_max_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "maximum", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_skew_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "skewness", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_kurt_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "kurtosis", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_median_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "median", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc05_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc05", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc25_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc25", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc75_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc75", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sig_perc95_xy"] = \
        lambda m : read_simple_table_by_key ("sigma_xx", "perc95", 6,
                        subtable_by_key(tables["Stress"],"m",str(m)))
lookups["sed_avg"] = \
        lambda : get_labelled_value ("average",
                                        tables["Strain Energy Density"])
lookups["sed_min"] = \
        lambda : get_labelled_value ("minimum",
                                        tables["Strain Energy Density"])
lookups["sed_max"] = \
        lambda : get_labelled_value ("maximum",
                                        tables["Strain Energy Density"])
lookups["sed_stddev"] = \
        lambda : get_labelled_value ("std_dev",
                                        tables["Strain Energy Density"])
lookups["sed_skew"] = \
        lambda : get_labelled_value ("skewness",
                                        tables["Strain Energy Density"])
lookups["sed_kurt"] = \
        lambda : get_labelled_value ("kurtosis",
                                        tables["Strain Energy Density"])
lookups["sed_median"] = \
        lambda : get_labelled_value ("median",
                                        tables["Strain Energy Density"])
lookups["sed_perc05"] = \
        lambda : get_labelled_value ("perc05",
                                        tables["Strain Energy Density"])
lookups["sed_perc25"] = \
        lambda : get_labelled_value ("perc25",
                                        tables["Strain Energy Density"])
lookups["sed_perc75"] = \
        lambda : get_labelled_value ("perc75",
                                        tables["Strain Energy Density"])
lookups["sed_perc95"] = \
        lambda : get_labelled_value ("perc95",
                                        tables["Strain Energy Density"])
lookups["sed_avg_mat%m"] = \
        lambda m : get_labelled_value ("average",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_min_mat%m"] = \
        lambda m : get_labelled_value ("minimum",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_max_mat%m"] = \
        lambda m : get_labelled_value ("maximum",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_stddev_mat%m"] = \
        lambda m : get_labelled_value ("std_dev",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_skew_mat%m"] = \
        lambda m : get_labelled_value ("skewness",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_kurt_mat%m"] = \
        lambda m : get_labelled_value ("kurtosis",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_median_mat%m"] = \
        lambda m : get_labelled_value ("median",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_perc05_mat%m"] = \
        lambda m : get_labelled_value ("perc05",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_perc25_mat%m"] = \
        lambda m : get_labelled_value ("perc25",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_perc75_mat%m"] = \
        lambda m : get_labelled_value ("perc75",
                        subtables(tables["Strain Energy Density"])[m])
lookups["sed_perc95_mat%m"] = \
        lambda m : get_labelled_value ("perc95",
                        subtables(tables["Strain Energy Density"])[m])
lookups["svm_avg"] = \
        lambda : get_labelled_value ("average",
                                        tables["Von Mises Stress"])
lookups["svm_min"] = \
        lambda : get_labelled_value ("minimum",
                                        tables["Von Mises Stress"])
lookups["svm_max"] = \
        lambda : get_labelled_value ("maximum",
                                        tables["Von Mises Stress"])
lookups["svm_stddev"] = \
        lambda : get_labelled_value ("std_dev",
                                        tables["Von Mises Stress"])
lookups["svm_skew"] = \
        lambda : get_labelled_value ("skewness",
                                        tables["Von Mises Stress"])
lookups["svm_kurt"] = \
        lambda : get_labelled_value ("kurtosis",
                                        tables["Von Mises Stress"])
lookups["svm_median"] = \
        lambda : get_labelled_value ("median",
                                        tables["Von Mises Stress"])
lookups["svm_perc05"] = \
        lambda : get_labelled_value ("perc05",
                                        tables["Von Mises Stress"])
lookups["svm_perc25"] = \
        lambda : get_labelled_value ("perc25",
                                        tables["Von Mises Stress"])
lookups["svm_perc75"] = \
        lambda : get_labelled_value ("perc75",
                                        tables["Von Mises Stress"])
lookups["svm_perc95"] = \
        lambda : get_labelled_value ("perc95",
                                        tables["Von Mises Stress"])
lookups["svm_avg_mat%m"] = \
        lambda m : get_labelled_value ("average",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_min_mat%m"] = \
        lambda m : get_labelled_value ("minimum",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_max_mat%m"] = \
        lambda m : get_labelled_value ("maximum",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_stddev_mat%m"] = \
        lambda m : get_labelled_value ("std_dev",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_skew_mat%m"] = \
        lambda m : get_labelled_value ("skewness",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_kurt_mat%m"] = \
        lambda m : get_labelled_value ("kurtosis",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_median_mat%m"] = \
        lambda m : get_labelled_value ("median",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_perc05_mat%m"] = \
        lambda m : get_labelled_value ("perc05",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_perc25_mat%m"] = \
        lambda m : get_labelled_value ("perc25",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_perc75_mat%m"] = \
        lambda m : get_labelled_value ("perc75",
                        subtables(tables["Von Mises Stress"])[m])
lookups["svm_perc95_mat%m"] = \
        lambda m : get_labelled_value ("perc95",
                        subtables(tables["Von Mises Stress"])[m])
lookups["dx_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "average", 1,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dx_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "std_dev", 1,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dx_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "minimum", 1,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dx_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "maximum", 1,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dx_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "median", 1,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dy_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "average", 2,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dy_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "std_dev", 2,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dy_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "minimum", 2,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dy_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "maximum", 2,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dy_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "median", 2,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dz_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "average", 3,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dz_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "std_dev", 3,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dz_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "minimum", 3,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dz_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "maximum", 3,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["dz_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("ux", "median", 3,
                        subtable_by_key(tables["Nodal Displacements"],"Node set",str(n)))
lookups["fx_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "total", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "average", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "std_dev", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "minimum", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "maximum", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "median", 1,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "total", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "average", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "std_dev", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "minimum", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "maximum", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fy_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "median", 2,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "total", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_avg_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "average", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_stddev_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "std_dev", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_min_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "minimum", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_max_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "maximum", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fz_median_ns%n"] = \
        lambda n : read_simple_table_by_key ("Fx", "median", 3,
                        subtable_by_key(tables["Nodal Forces"],"Node set",str(n)))
lookups["fx_ns%n_mat%m"] = \
            lambda n,m : read_simple_table_by_row ("material", m, 1,
                        subtable_by_key(tables["Load Sharing"],"Node set",str(n)))
lookups["fy_ns%n_mat%m"] = \
            lambda n,m : read_simple_table_by_row ("material", m, 2,
                        subtable_by_key(tables["Load Sharing"],"Node set",str(n)))
lookups["fz_ns%n_mat%m"] = \
            lambda n,m : read_simple_table_by_row ("material", m, 3,
                        subtable_by_key(tables["Load Sharing"],"Node set",str(n)))
lookups["pis_stiffx"] = \
            lambda : get_labelled_triplet ("Axial stiffness",
                                            tables["Pistoia Failure Load Estimate"])[0]
lookups["pis_stiffy"] = \
            lambda : get_labelled_triplet ("Axial stiffness",
                                            tables["Pistoia Failure Load Estimate"])[1]
lookups["pis_stiffz"] = \
        lambda : alternates ( 
            [ lambda : get_labelled_triplet ("Axial stiffness",
                                            tables["Pistoia Failure Load Estimate"])[2],
            lambda : get_labelled_value ("Axial stiffness",
                                            tables["Pistoia Failure Load Estimate"]),
            lambda : get_labelled_value ("Axial stiffness \[N/mm\]",
                                            tables["Pistoia Failure Load Estimate"]) ])
lookups["pis_fx_fail"] = \
            lambda : get_labelled_triplet ("Failure load \(RF \* factor\)",
                                            tables["Pistoia Failure Load Estimate"])[0]
lookups["pis_fy_fail"] = \
            lambda : get_labelled_triplet ("Failure load \(RF \* factor\)",
                                            tables["Pistoia Failure Load Estimate"])[1]
lookups["pis_fz_fail"] = \
        lambda : alternates ( 
            [ lambda : get_labelled_triplet ("Failure load \(RF \* factor\)",
                                            tables["Pistoia Failure Load Estimate"])[2],
            lambda : get_labelled_value ("Failure load \(RFz \* factor\)",
                                            tables["Pistoia Failure Load Estimate"]),
            lambda : get_labelled_value ("Failure load \(RFz \* factor\) \[N\]",
                                            tables["Pistoia Failure Load Estimate"]) ])
