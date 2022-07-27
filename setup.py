import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='pace',
        version="0.1",
        description='A Workflows approach to biochemical supramolecular polymorphism',
        author='Srinivas Mushnoori',
        url='',
        license='',
        packages=setuptools.find_packages(),
        install_requires=['numpy',
                          'gitpython',
                          'radical.entk'],
        scripts=['src/pace2.py',
                 'src/candidate.py',
                 'src/candidate_manager.py'],
        zip_safe=True,
        )
