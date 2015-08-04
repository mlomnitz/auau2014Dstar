#!/bin/csh

#   Submit fileLists for classes derived from 
#    StPicoHFMaker;
#
#  - script will create a folder ${baseFolder}/jobs
#    all submission related files will end up there
#
#  - in ${baseFolder} the script expects (links or the actual folders)
#      .sl64_gcc447
#      StRoot                     ( from the git repo )
#      run14AuAu200GeVPrescales   ( from the git repo )
#      starSubmit                 ( from the git repo )
#
#      picoLists                  ( from the fileList git repo )
#
#   - the rootMacro is expected in StRoot/macros
#
#   - the bad run list is expected in ${baseFolder}
#     or in ${baseFolder}/picoLists
#
# ###############################################

# -- baseFolder of job
set baseFolder=/global/homes/m/mlomnitz/mlomnitz_projectdir/Dstar

# --input file 
#    makerMode 0,1 : list must contain picoDst.root files
#    makerMode 2   : list must contain picoHFtree.root files
set input=${baseFolder}/RunQvector.list

# -- set maker mode
#    0 - kAnalyze, 
#    1 - kWrite
#    2 - kRead
set makerMode=0

# -- set root macro
set rootMacro=runPicoDstarAnaMaker.C

# -- set filename for bad run list
set badRunListFileName="picoList_bad_MB.list"

# ###############################################
# -- DON'T CHANGE BELOW THAT LINE
# ###############################################

# -- production Id (kAnalyse / kRead)
set productionId=`date +%F_%H-%M`

# -- production Id (kWrite)
set treeName=hfTree

# -- job submission directory
mkdir -p ${baseFolder}/jobs

# -- result directory
mkdir -p ${baseFolder}/production

pushd ${baseFolder}/jobs > /dev/null

# -- prepare folder
mkdir -p report err log list csh

# -- check for prerequisits and create links
set folders=".sl64_gcc447 StRoot run14AuAu200GeVPrescales starSubmit"

foreach folder ( $folders ) 
    if ( ! -d ${baseFolder}/${folder} ) then
	echo "${folder} does not exist in ${baseFolder}"
	exit
    else
	ln -sf  ${baseFolder}/${folder}
    endif
end

if  ( ! -e ${baseFolder}/StRoot/macros/${rootMacro} ) then
    echo "${rootMacro} does not exist in ${baseFolder}/StRoot/macros"
    exit
endif

if ( ! -e ${input} ) then
    echo "Filelist ${input} does not exist"
    exit
endif

if ( ! -e ${baseFolder}/starSubmit/submitPicoDstarAnaMaker.xml ) then
    echo "XML submitPicoDstarAnaMaker.xml does not exist"
    exit
else
    ln -sf ${baseFolder}/starSubmit/submitPicoDstarAnaMaker.xml 
endif

if  ( -e ${baseFolder}/${badRunListFileName} ) then
    cp  ${baseFolder}/${badRunListFileName} picoList_badRuns.list
else if ( -e ${baseFolder}/picoLists/${badRunListFileName} ) then
    cp  ${baseFolder}/picoLists/${badRunListFileName} picoList_badRuns.list
else
    echo "${badRunListFileName} does not exist in ${baseFolder} nor ${baseFolder}/picoLists"
    exit
endif

if ( -e LocalLibraries.zip ) then
    rm LocalLibraries.zip
endif 

if ( -d LocalLibraries.package ) then
    rm -rf LocalLibraries.package
endif 

# -- submit 
star-submit-template -template submitPicoDstarAnaMaker.xml -entities listOfFiles=${input},basePath=${baseFolder},prodId=${productionId},mMode=${makerMode},treeName=${treeName},rootMacro=${rootMacro}

popd > /dev/null
