FROM ubuntu:23.04

RUN apt update
RUN apt install -y vim curl wget g++ gcc cmake locate man libomp-dev zip unzip software-properties-common gcc-multilib bash python3
WORKDIR /usr/share/
RUN apt install -y libzmq5-dev git git-lfs libstdc++6
RUN git clone "https://github.com/arnabpune/SPOTLIGHT.git"
#RUN wget "https://onedrive.live.com/download?resid=736A06E59DB7B08E%213273&authkey=!ACVxFEnlsAhQisk" -O libtorch.zip
RUN wget "localhost:8000/libtorch.zip" -O libtorch.zip
RUN unzip libtorch.zip
WORKDIR /usr/share/SPOTLIGHT/spotlight_pt_port/
RUN ls -ltr --color=auto
RUN sed -e '0,/LIBTORCH_LOC/ {/LIBTORCH_LOC/d}' autobuild.sh > ab
RUN echo "LIBTORCH_LOC=/usr/share/libtorch" > autobuild.sh
RUN cat ab >> autobuild.sh
RUN mkdir -p /usr/share/SPOTLIGHT/models/default
RUN cp ../spotlight_data/model/*.pt /usr/share/SPOTLIGHT/models/default
RUN chmod +x setpaths.sh && ./setpaths.sh /usr/share/SPOTLIGHT/models/default
RUN chmod +x autobuild.sh && ./autobuild.sh
WORKDIR /usr/share/SPOTLIGHT/spotlight_pt_port/build
RUN ls && cp *PT* /usr/bin
RUN cp -n *production* /usr/bin
WORKDIR /usr/share/SPOTLIGHT/scripts
RUN chmod +x ./* && cp ./* /usr/bin

ENTRYPOINT ["/bin/bash","-l","-c"]
CMD ["/bin/bash"]
