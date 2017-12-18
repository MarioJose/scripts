#!/bin/bash
#
# Version: 1.0.0
#

###########
# VARIABLES
###########

# Regular expression to search date and version format
DATE_REGEX='[0-9]{4}\-[0-9]{2}\-[0-9]{2}'
VERSION_REGEX='[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}'

# Path to package
PKG_PATH='./my_package/'

# A list of files to update
FILES2UPDATE=('DESCRIPTION' 'README.md' 'NEWS' 'man/my_package-package.Rd')

# Default values of variables
DATE_NEW=''
VERSION_NEW=''
DATE_OLD=`cat ${PKG_PATH}DESCRIPTION | grep -oE $DATE_REGEX`
VERSION_OLD=`cat ${PKG_PATH}DESCRIPTION | grep -oE $VERSION_REGEX`
## Flag to add version date in NEWS file
ADD_IN_NEWS=0


###########
# FUNCTIONS
###########

function checkVersion {
    if [ "$1" != "" ]
    then
        CHK_VERSION=`echo $1 | grep -oE $VERSION_REGEX`
        if [ "$CHK_VERSION" != "" ]
        then
            VERSION_NEW=$1
        else
            printf 'Version do not have format NUMBER.NUMBER.NUMBER \n'
            exit 1
        fi
    else
        printf 'Version not specified\n'
        exit 1
    fi

    if [ "$VERSION_NEW" == "$VERSION_OLD" ]
    then
        printf '\nNew version is equal to old version \n'
        read -p 'Do you want continue (y/n)?' USER_CHOICE
        case "$USER_CHOICE" in
            y|Y ) ;;
            n|N ) exit 1;;
            * );;
        esac
    fi
}

function checkDate {
    if [ "$1" != "" ]
    then
        CHK_DATE=`echo $1| grep -oE $DATE_REGEX`
        if [ "$CHK_DATE" != "" ]
        then
            DATE_NEW=$1
        else
            printf 'Date do not have format YYYY-MM-DD \n'
            exit 1
        fi
    else
        DATE_NEW=`date +%Y-%m-%d`
    fi

    if [ "$DATE_NEW" == "$DATE_OLD" ]
    then
        printf '\nNew date is equal to old date \n'
        read -p 'Do you want continue (y/n)?' USER_CHOICE
        case "$USER_CHOICE" in
            y|Y ) ;;
            n|N ) exit 1;;
            * );;
        esac
    fi
}

function showInfo {
    for i in ${FILES2UPDATE[@]}
    do
        printf "\nFile: ${i}\n"
        echo "Date:" `cat ${PKG_PATH}${i} | grep -m 1 -oE $DATE_REGEX`
        echo "Version:" `cat ${PKG_PATH}${i} | grep -m 1 -oE $VERSION_REGEX`
    done
}

function updateInfo {
    for i in ${FILES2UPDATE[@]}
    do
        printf "\nFile: ${i}\n"
        echo "Updating date from $DATE_OLD to $DATE_NEW"
        sed -i -E "0,/$DATE_REGEX/ s//$DATE_NEW/" ${PKG_PATH}${i}
        echo "Updating version from $VERSION_OLD to $VERSION_NEW"
        sed -i -E "0,/$VERSION_REGEX/ s//$VERSION_NEW/" ${PKG_PATH}${i}
     done
}

function showUsage {
    printf "Usage: $0 -v version (required) -d date (optional) -n (optional)\n"
    printf " -v: New version in N.N.N format (N is numerical);\n"
    printf " -d: New date in YYYY.MM.DD format. If not informed,\n     it is set to current date;\n"
    printf " -n: New date and version entry in NEWS file.\n     If -n is not informed, date and version will be updated in NEWS.\n"
}

function addInNews {
    printf "\nAdding version and date in NEWS file\n"
    sed -i -E "1s/^/$VERSION_NEW \($DATE_NEW\)\n\*\ \n\n/" "${PKG_PATH}NEWS"
}

function parseArguments {
    while [[ $# -gt 0 ]]
    do
        case $1 in
            -v ) checkVersion $2
                 shift 2
                 ;;
            -d ) checkDate $2
                 shift 2
                 ;;
            -n ) ADD_IN_NEWS=1
                 shift
                 ;;
            * ) showUsage
                exit 1
        esac
    done

    if [ "${VERSION_NEW}" == "" ]
    then
        showUsage
        exit 1
    fi
    
    if [ "${DATE_NEW}" == "" ]
    then
        checkDate
    fi
}


######
# MAIN
######

if [ $# -eq 0 ]
then
    showUsage
    exit 1
else
    parseArguments "$@"
fi

if [ $ADD_IN_NEWS -eq 1 ]
then
    for i in "${!FILES2UPDATE[@]}"
    do
        if [ "${FILES2UPDATE[$i]}" == "NEWS" ]
        then
            unset 'FILES2UPDATE[$i]'
        fi
    done
fi

showInfo
printf "\nNew values\n"
echo "Version: $VERSION_NEW"
echo "Date: $DATE_NEW"
printf "\n"
read -p 'Update to new values (y/n)?' USER_CHOICE
case "$USER_CHOICE" in
    y|Y ) updateInfo ;;
    n|N ) ;;
    * ) exit 0;;
esac

if [ $ADD_IN_NEWS -eq 1 ]
then
    printf "\nFile: NEWS (first 10 lines)\n"
    sed -n 1,10p "${PKG_PATH}NEWS"
    
    CHK_NEWS=`grep -oE "${VERSION_NEW}\ \(${DATE_NEW}\)" ${PKG_PATH}NEWS`
    if [ "${CHK_NEWS}" == "" ]
    then
        printf "\n"
        read -p 'Add date-version to NEWS file (y/n)?' USER_CHOICE
        case "$USER_CHOICE" in
            y|Y ) addInNews ;;
            n|N ) ;;
            * ) exit 0;;
        esac
    else
        printf "\nDate-version already in NEWS file\n"
    fi
fi




