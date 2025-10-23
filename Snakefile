#snakeobjects 3.1.2
#Snakeobjects types: ['chunk', 'in_chromosome']
shell.prefix("set -euo pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, EF, add_targets, get_targets
from snakeobjects.remoteProjects import download_project_files_from_remote

# CLOUD related reorganization
import os,sys

if ("SO_CONTAINER" in os.environ and 
  "SO_REMOTE" in os.environ and
  os.environ["SO_CONTAINER"] == 'yes'):
    if ":" in os.environ["SO_REMOTE"]:
      provider, bucket = os.environ["SO_REMOTE"].split(":")
      download_project_files_from_remote(provider,bucket)
      os.environ['SO_PIPELINE'] = "workflow"
      os.environ['SO_PROJECT'] = "./" + bucket 
    fixPP = os.environ['SO_PIPELINE'] + "/fix_permissions.sh"
    if os.path.exists(fixPP):
      os.system("bash %s" % fixPP)

project = Project()
so_pipeline=project.get_pipeline_directory()

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

include: "chunk.snakefile"
rule so_all_chunk:
  input:
    expand("{of}", of=project.get_all_object_flags("chunk"))

rule so_chunk_obj:
  input: get_targets("chunk")
  output: touch("chunk/{oid}/obj.flag")

include: "in_chromosome.snakefile"
rule so_all_in_chromosome:
  input:
    expand("{of}", of=project.get_all_object_flags("in_chromosome"))

rule so_in_chromosome_obj:
  input: get_targets("in_chromosome")
  output: touch("in_chromosome/{oid}/obj.flag")

