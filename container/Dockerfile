FROM dolfinx/lab:stable
ENV TZ=Asia/Tokyo
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata
RUN apt-get install -y dnsutils
RUN apt-get install -y net-tools
RUN apt-get install -y x11-apps
RUN apt-get install -y vim nano
RUN apt-get install -y lsof
RUN apt-get install -y xdg-utils
RUN apt-get install -y netcat
RUN apt-get install -y gtkwave
RUN mkdir /usr/share/desktop-directories
RUN pip install pandas
RUN pip install h5py
RUN pip install 'jupyterlab>=3' ipywidgets 'pyvista[all]' ipyvtklink
RUN pip install pcb-tools
RUN pip install cairocffi
RUN pip install gmsh
