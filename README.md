# Computational_docker_protocol

## Docker Images
- __md_protocol__  [![Docker Pulls](https://img.shields.io/docker/pulls/ychen209/md_protocol)](https://hub.docker.com/r/ychen209/md_protocol)  
An Image used for molecular dynamics simulation using OpenMM, including protein-based and protein-ligand complex.

- __md_ramachandran_protocol__  [![Docker Pulls](https://img.shields.io/docker/pulls/ychen209/md_ramachandran)](https://hub.docker.com/r/ychen209/md_ramachandran)  
An Image used for Ramachandran plot of $\alpha$-amino acids (natural and unnatural), in which the phi/psi angles were calculated by molecular dynamics simulation.
# run time example (default 100 ns per uaa): RTX 4090 ~ 1300 ns, RTX 4060 laptop ~ 530 ns, RTX 3090 ~ 1100 ns, RTX 3080 Ti ~ 1100 ns, Tesla T4 ~ 660 ns.

## Fundamental Docker Usage

### Key Concepts  
- __Image__: Application + environment template  
- __Container__: Running instance of an image  
- __Dockerfile__: Instructions to build an image  
- __Registry__: Stores and distributes images  

### Installation  
Following the __Manuals__ from [Docker](https://www.docker.com/).  
Docker supports a wide range of platforms across operating systems (including Windows, macOS and Linux), architectures, and cloud environments.

### Common Commands  
The common used commands are listed as following:  
```
# Image
docker pull <image>           # Pull an image
docker images                 # List local images
docker rmi <image_id>         # Remove an image

# Container
docker run <image>            # Run an image
docker ps                     # List running containers
docker ps -a                  # List all containers
docker stop <container_id>    # Stop a container
docker rm <container_id>      # Remove a container

# Dockerfile
FROM python:3.10-slim
WORKDIR /app
COPY . .
RUN pip install -r requirements.txt
CMD ["python", "app.py"]

# Build from Dockerfile
docker build -t <myImage> .
```

### Some options for `docker run`  
| Option | Description |
|--------|-------------|
| `--rm` | Automatically remove the container after it exits. |
| `-i` | Keep STDIN open. |
| `-t` | Allocate a pseudo-TTY (use with `-i` → `-it`). |
| `--env-file <file>` | Load environment variables from a file. |
| `-v <host:container>` | Bind mount a host directory or file. |
| `--mount` | Advanced mount configuration. |
| `-m`, `--memory <limit>` | Memory limit for the container. |
| `--cpus <count>` | Limit CPU usage. |
| `--cpu-shares <value>` | Set CPU scheduling weight. |
| `--gpus <all/N>` | Assign GPU devices. |
| `--user <uid:gid>` | Run container as a specific user. |
| `--privileged` | Give the container full privileges. |
| `--pull <policy>` | Pull image policy (`always`, `missing`, `never`). |

&nbsp;

## md_ramachandran_protocol

### Input  
The input file should contain two columns with header of `name` and `smiles`. The file should be in the format of comma-separated values.  
The following is an example:
```
# example file name as 'aa_file.csv'
name,smiles
acK,CC(=O)NCCCC[C@@H](C(=O)O)N
bocK,CC(C)(C)OC(=O)C(CCC[C@@H](C(=O)O)N)N
```

### Run  
```
# example working directory
workdir/
└── input/
      └── aa_file.csv             # input uaa file (only needed file for input)
```

```
# run the command in your working directory for example 'workdir'
# usage
docker run -it --rm --gpus all -v {local input folder}:{container input folder} -v {local output folder}:{container output folder} ychen209/md_ramachandran --uaa_file {uaa_file} --output {output folder}

{local input folder} - the input folder in your working directory, './input' here
{container input folder} - the input folder path in the container, and Docker mounts the {local input folder} to this internal path, '/input' recommended here
{local output folder} - the output folder in your working directory and the folder will be create if not exist, './output' recommended here
{container output folder} - the output folder path in the container, and Docker mounts the {local output folder} to this internal path, '/output' recommended here
{uaa_file} - the file containing uaa smiles for calculation in your working directory, 'input/aa_file.csv' here
{output folder} - the folder for output, '/output' recommended here

# default example
docker run -it --rm --gpus all -v ./input:/input -v ./output:/output ychen209/md_ramachandran --uaa_file /input/aa_file.csv --output /output

# if change the parameters
docker run -it --rm --gpus all -v ./input:/input -v ./output:/output ychen209/md_ramachandran --uaa_file /input/aa_file.csv --output /output --timestep 0.002 --ministep 50000 --steps 500000 --savesteps 500
```
### Options
| Option <default> | Description |
|---|---|
| `--uaa_file <>` | The file containing uaa `name` and `smiles` for calculation in the `csv` format. |
| `--output <.>` | Assign an output folder. The folder will be created if not exist. |
| `--timestep <0.002>` | Assign the time duration of each step in the `* unit of picosecond` format. This parameter is not recommended to modify, while increase of this parameter will highly increase computational load and evenleading to crash. |
| `--ministep <500000>` | Steps for minimization. Default is 1 ns. |
| `--steps <50000000>` | Steps for simulation. Default is 100 ns. |
| `--savesteps <50000>` | Steps for saving each frame. Default is 100 ps per frame. |

### Example output
```
# acK as example
workdir/
├── input/
│     └── aa_file.csv             # input uaa file (only needed file for input)
└── output/
      ├── log.txt                 # status of the batch
      ├── error.txt               # generated if some uaa failed
      ├── md/
      │    ├── acK.sdf            # generated tripeptide
      │    ├── acK.cif            # topology structure in md   
      │    ├── acK.dcd            # trajectory file in md
      │    └── acK_log.txt        # information during md
      └── rama/
           ├── acK.txt            # phi/psi information during md
           ├── acK.png            # Ramachandran plots of the free energy surface
           ├── acK.svg            # plot in vector graphics
           ├── acK_scatter.png    # plot with scatter
           └── acK_scatter.svg    # plot with scatter in vector graphics
```
