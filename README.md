# MCP

Microbiome Coverage Profiler

In the repository there is a virtual environment (venv) with python and the necessary packages installed already
Access the virtual environment using the command

`source venv/bin/activate`

mbprofile.py takes in 2 arguments:
  1. The csv file with mapped reads (str)
  2. the read length (int)

To run the code you can call it like this:

`python mbprofile.py rop_output_file.csv 100`

After you are done with the virtual environment, you can close the environment with the command

`deactivate`
