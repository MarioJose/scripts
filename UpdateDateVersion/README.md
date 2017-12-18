# Bash script to update date and version in R package

### Version: 1.0.0

`update_date_version.sh` script check date and version of an release in informed files. It will check format and update date and version of files in array variable `FILES2UPDATE`.

In `update_date_version.sh` you must define `PKG_PATH` and `FILES2UPDATE` variables.

    * `PKG_PATH`: String with path to package folder;
    
    * `FILES2UPDATE`: Array with files to update.

By default, version must be N.N.N format, with _N_ as number, and date with YYYY-MM-DD format. This is specified as regular expression in `DATA_REGEX` and `VERSION_REGEX` variables.

## Usage
Usage: ./update\_date\_version.sh -v version (required) -d date (optional) -n (optional)

 -v: New version in N.N.N format (N is numerical);
 
 -d: New date in YYYY.MM.DD format. If not informed, it is set to current date;
 
 -n: New date and version entry in NEWS file. If -n is not informed, date and version will be updated in NEWS.


