lsetup "views LCG_97python3 x86_64-centos7-gcc8-opt"

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
PYTHONPATH=$PYTHONPATH:$parent_path/core_code

LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF
