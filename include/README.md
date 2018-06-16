This folder is segmented as follows:
  - `latmrg/` contains the headers describbing the low-level API of LatMRG.
  - `latmrg/mrgtypes` contains a few specific implementations of MRG or MMRG 
    that can be used without reimplementation.
  - `latmrg-high` contains high level routines of the LatMRG API.

This segmentation has two goals.
  - First, it makes the directories easier to navigate and let us dodge the
    difficulty of managing/finding a file in a 50+ files directory
  - Second, it makes the API more accessible. An unexperienced user will know 
    that he has to first look at the routines availaible in the modules of
    `latmrg-high` because those are meant to be used externally.
