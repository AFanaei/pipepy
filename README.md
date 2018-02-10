# pipepy #
[![Build Status](https://travis-ci.org/jajool/pipepy.svg?branch=master)](https://travis-ci.org/jajool/pipepy)


a library for simulating pipeline.
## how to run project: ##
VScode:
* install VScode.
* install python extention for vscode.
* update file > preferences > settings:
```
{
    "python.linting.pylintEnabled": false,

    "python.linting.flake8Enabled": true,

    "files.exclude": {
        "**/*.pyc": true
    }
}
```

## git ##
* install git
* register in https://github.com/
* create new ssh key and add to github https://help.github.com/articles/connecting-to-github-with-ssh/
* git clone <REPOSITORY ADDRESS> # get files 
* git commit # save changes.
* git push origin master # send changes to server.
* git pull origin master # for getting latest changes.


## conda ##
* anaconda commands: https://conda.io/docs/_downloads/conda-cheatsheet.pdf
* conda env create --name pipepy --file bio-env.txt 
* activate pipepy
* python test.py
