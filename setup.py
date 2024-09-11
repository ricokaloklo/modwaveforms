import setuptools

verstr = "0.0.1"

setuptools.setup(
    name="modwaveforms",
    version=verstr,
    author="Jose María Ezquiaga",
    author_email="jose.ezquiaga@nbi.ku.dk",
    description="Collection of modified waveforms for parameter estimation",
    packages=[
        "modwaveforms",
    ],

     install_requires=[
        "bilby",
        "bilby_pipe",
    ],

    classifiers=[
        "Programming Language :: Python :: 3.7",
    ],
    python_requires='>=3.7',
)
