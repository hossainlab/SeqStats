import setuptools

# Open README File 
with open("README.md", "r") as fh:
    long_description = fh.read()

# Initial Setup 
setuptools.setup(
    # Module / Package Name 
    name="SeqStats", 
    # Version of Package 
    version="0.0.1",
    # Name of Author 
    author="Jubayer Hossain",
    # Author Email 
    author_email="contact.jubayerhossain@gmail.com",
    # Small Description about Package / Module 
    description="SeqStats-A Bioinformatics Package for Sequence Analysis",
    # Long Description of Package 
    long_description="SeqStats - A Bioinformatics Package for DNA Sequence Analysis",
    # Specifying taht we are using markdown file for description 
    long_description_content_type="text/markdown",
    # Any link to reach this module, if you have any webpage or github profile 
    url="https://github.com/hossainlab/pydemo",
    packages=setuptools.find_packages(),

    # Classifiers like program is suitable for python3, just leave as it is. 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

