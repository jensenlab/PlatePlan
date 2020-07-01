import setuptools

setuptools.setup(
    name="plateplan",
    version="0.0.1",
    author="Adam Dama, Paul Jensen",
    author_email="software@jensenlab.net",
    description="PlatePlan is meant to be the bridge between your hypotheses and results. In order to make that jump, a lot goes into the logistics and planning out an experiment, from how many multi-well plates to use and where to dispense each reagent, to which reagents to batch together and how much volume to dispense. This logistics are especially challenging in the case of high-throughput assays where you may have tens of thousands of experiments to manage and keep track of.",
    packages=setuptools.find_packages(),
    install_requires=["numpy", "pandas", "pyparsing"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
