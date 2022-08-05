import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
     name='hpsandbox',  
     version='2.0.0',
     scripts=['hpsandbox'] ,
     author="Vincent Voelz",
     author_email="vvoelz@gmail.com",
     description="Python package to explore the 2D HP model of Chan and Dill.",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/vvoelz/HPSandbox",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )
