from morphomics import protocols
import tomli
import sys

Parameters_ID = int(sys.argv[1])
parameters_filepath = (
    "/media/drive_siegert/RyCu/Projects/Git_codes/morphomics_v2/Morphomics.Parameters.%d.toml"%Parameters_ID
)
with open(parameters_filepath, mode="rb") as _parameter_file:
    parameters = tomli.load(_parameter_file)
parameters["Parameters_ID"] = Parameters_ID

protocol = protocols.Protocols(parameters, Parameters_ID)
script_sequence = parameters["Protocols"]
for sequence in script_sequence:
    print("Doing %s..."%sequence)
    perform_this = getattr(protocol, sequence)
    perform_this()