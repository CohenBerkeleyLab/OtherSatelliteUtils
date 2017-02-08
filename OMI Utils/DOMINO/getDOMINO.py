from __future__ import print_function
import subprocess
import urllib
import datetime as dt
import os
import sys
import argparse
import re

root_url = "http://www.temis.nl/airpollution/no2col/data/omi/data_v2/%Y/"
file_pattern = "omi_no2_he5_%Y%m%d.tar"

def __shell_error(msg, exitcode=1):
    print(msg, file=sys.stderr)
    exit(exitcode)

def get_args():
    parser = argparse.ArgumentParser(description="Downloads OMI DOMINO files to the current directory for the given date range")
    parser.add_argument("startdate", help="First day to download in yyyy-mm-dd format")
    parser.add_argument("enddate", help="Last day to download in yyyy-mm-dd format")

    args = parser.parse_args()
    date_restr = "\d\d\d\d-\d\d-\d\d"
    if not re.match(date_restr, args.startdate):
        __shell_error("startdate is not in yyyy-mm-dd format")
    elif not re.match(date_restr, args.enddate):
        __shell_error("enddate is not in yyyy-mm-dd format")

    sdate = dt.datetime.strptime(args.startdate, "%Y-%m-%d")
    edate = dt.datetime.strptime(args.enddate, "%Y-%m-%d")

    if edate < sdate:
        __shell_error("startdate cannot be later than enddate")

    return sdate, edate

def download_and_unpack_file(filedate):
    # Download the file
    year_str = "{:04}".format(filedate.year)
    month_str = "{:02}".format(filedate.month)
    day_str = "{:02}".format(filedate.day)
    this_url = root_url.replace("%Y", year_str)
    this_file = file_pattern.replace("%Y", year_str).replace("%m", month_str).replace("%d", day_str)
    print("Downloading ", this_url + this_file)
    result = urllib.urlretrieve(this_url + this_file, this_file)

    # Untar the archive
    subprocess.check_output(["tar", "-xf", this_file])
    os.remove(this_file)

def main():
    sdate, edate = get_args()
    td = dt.timedelta(days=1)
    curr_date = sdate
    while curr_date <= edate:
        download_and_unpack_file(curr_date)
        curr_date += td

if __name__ == "__main__":
    main()