import pathlib

import setuptools

__packagename__ = 'catibts'
ROOT = pathlib.Path(__file__).parent


def get_version():
    import os, re
    VERSIONFILE = os.path.join(__packagename__, '__init__.py')
    initfile_lines = open(VERSIONFILE, 'rt').readlines()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    for line in initfile_lines:
        mo = re.search(VSRE, line, re.M)
        if mo:
            return mo.group(1)
    raise RuntimeError('Unable to find version string in %s.' % (VERSIONFILE,))

__version__ = get_version()


setuptools.setup(name = __packagename__,
        packages = setuptools.find_packages(), # this must be the same as the name above
        version=__version__,
        description='Package for Thompson Sampling of CatIB producers.',
        url='https://github.com/JuBiotech/Supplement-to-Helleckes-Kuesters-et-al.-2023',
        author='Laura Helleckes',
        copyright='(c) 2023 Forschungszentrum Jülich GmbH',
        license='(c) 2023 Forschungszentrum Jülich GmbH',
        classifiers= [
            'Programming Language :: Python',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.8',
            'Intended Audience :: Developers'
        ],
        install_requires=open(pathlib.Path(ROOT, "requirements.txt")).readlines()
)

