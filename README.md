# Computational_docker_protocol

## Docker Images
- __md_protocol__  [![Docker Pulls](https://img.shields.io/docker/pulls/ychen209/md_protocol)](https://hub.docker.com/r/ychen209/md_protocol)  
An Image used for molecular dynamics simulation using OpenMM, including protein-based and protein-ligand complex.

- __md_ramachandran_protocol__  [![Docker Pulls](https://img.shields.io/docker/pulls/ychen209/md_ramachandran)](https://hub.docker.com/r/ychen209/md_ramachandran)  
An Image used for Ramachandran plot of $\alpha$-amino acids (natural and unnatural), in which the phi/psi angles were calculated by molecular dynamics simulation.

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

### Options for `docker run`  
| Option | Description |
|--------|-------------|
| `-d`, `--detach` | Run container in the background. |
| `--rm` | Automatically remove the container after it exits. |
| `--name <name>` | Assign a custom container name. |
| `-i` | Keep STDIN open. |
| `-t` | Allocate a pseudo-TTY (use with `-i` → `-it`). |
| `-p <host:container>` | Map host port to container port. |
| `-P` | Publish all exposed ports automatically. |
| `--network <network>` | Connect container to a specific Docker network. |
| `--dns <server>` | Set DNS server for the container. |
| `-e KEY=value` | Set environment variable. |
| `--env-file <file>` | Load environment variables from a file. |
| `-v <host:container>` | Bind mount a host directory or file. |
| `--mount` | Advanced mount configuration. |
| `--read-only` | Set container filesystem as read-only. |
| `--tmpfs <path>` | Mount a tmpfs (memory-only filesystem). |
| `-m`, `--memory <limit>` | Memory limit for the container. |
| `--cpus <count>` | Limit CPU usage. |
| `--cpu-shares <value>` | Set CPU scheduling weight. |
| `--gpus <all/N>` | Assign GPU devices. |
| `--restart <policy>` | Set restart policy (`no`, `on-failure`, `always`, `unless-stopped`). |
| `--user <uid:gid>` | Run container as a specific user. |
| `--privileged` | Give the container full privileges. |
| `--cap-add <cap>` | Add Linux capabilities. |
| `--cap-drop <cap>` | Drop Linux capabilities. |
| `--security-opt <opt>` | Set security options (AppArmor, seccomp, etc.). |
| `--entrypoint <cmd>` | Override the image entrypoint. |
| `-w`, `--workdir <dir>` | Set working directory inside the container. |
| `--pull <policy>` | Pull image policy (`always`, `missing`, `never`). |

&nbsp;

## md_ramachandran_protocol

### Input  
The input file should contain two columns with header of `name` and `smiles`. The file should be in the format of comma-separated values.  
The following is an example:
```
name,smiles
acK,CC(=O)NCCCC[C@@H](C(=O)O)N
bocK,CC(C)(C)OC(=O)C(CCC[C@@H](C(=O)O)N)N
```

### Run  
```
# default
docker run -it --rm --gpus all -v ./input:/input -v ./output:/output ychen209/md_ramachandran --uaa_file /input/uaa.txt --output /output

# example
docker run -it --rm --gpus all -v ./input:/input -v ./output:/output ychen209/md_ramachandran --uaa_file /input/uaa.txt --output /output --timestep 0.002 --ministep 50000 --steps 500000 --savesteps 500
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

### Output
```
# acK as example
workdir/
├── input/
│     └── uaa.txt            # input uaa file
└── output/
      ├── md/
      │    ├── acK.sdf       # generated tripeptide
      │    ├── acK.cif       # topology structure in md   
      │    ├── acK.dcd       # trajectory file in md
      │    └── acK_log.txt   # information during md
      └── rama/
           ├── acK.txt       # phi/psi information during md
           ├── acK.png       # Ramachandran plot
           └── acK.svg       # Ramachandran plot in vector graphics
```
