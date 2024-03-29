from setuptools import setup
from os import path as os_path


## INFOS ##
package     = 'doebase'
descr       = 'Optimal Design of Experiments (DoE) base package for synthetic biology'
url         = 'https://github.com/pablocarb/doebase'
authors     = 'Pablo Carbonell'
corr_author = 'pablo.carbonell@upv.es'

## LONG DESCRIPTION
with open(
    os_path.join(
        os_path.dirname(os_path.realpath(__file__)),
        'README.md'
    ),
    'r',
    encoding='utf-8'
) as f:
    long_description = f.read()

def get_version():
    with open(
        os_path.join(
            os_path.dirname(os_path.realpath(__file__)),
            'CHANGELOG.md'
        ),
        'r'
    ) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('##'):
            from re import search
            m = search("\[(.+)\]", line)
            if m:
                return m.group(1)

setup(
    name                          = package,
    version                       = get_version(),
    author                        = authors,
    author_email                  = corr_author,
    description                   = descr,
    long_description              = long_description,
    long_description_content_type = 'text/markdown',
    url                           = url,
    packages                      = [package],
    package_dir                   = {package: package},
    include_package_data          = True,
    test_suite                    = 'pytest',
    license                       = 'MIT',
    classifiers                   = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires               = '>=3.7',
)

