# borrowed from cobrapy package https://github.com/opencobra/cobrapy/blob/devel/cobra/clean.sh
find . -type f -regex '.*pyc' | xargs rm
find . -type f -regex '.*class' | xargs rm
