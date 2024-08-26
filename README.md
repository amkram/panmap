## Installation
**Docker local**

`git clone https://github.com/amkram/panmap`

`cd panmap && docker build -t panmap .`

`docker run -it -v .:/panmap panmap bash`

The `-v` argument above mounts the current directory (panmap source) to `/panmap` in the Docker container (`/panmap` is the source directory within the container). This means that changes in your local panmap source dir are reflected to the container's panmap source dir. You could also separate them, e.g., `-v .:/mine`, which links `.` to `/mine` in the container, leaving `/panmap` alone.
