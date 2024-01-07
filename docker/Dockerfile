FROM ubuntu:23.04
ENTRYPOINT ["/bin/bash"]

RUN apt update
RUN apt install -y vim curl wget g++ gcc cmake locate man libomp-dev zip unzip software-properties-common gcc-multilib
WORKDIR /usr/share/
RUN apt install -y libzmq5-dev git git-lfs libstdc++6
RUN git clone "https://github.com/arnabpune/SPOTLIGHT.git"
RUN wget "https://public.bl.files.1drv.com/y4m5BgKphUayNZHbMeatJ2p8R0xShxoUKgV5SGJojqFuLY83jRcq_1-nq1OhHLIpRH2kg9y-_RADZbqyeuden68Uidps3mTZvRnWTTaPM24ASasJoRYARYXS1WPjxazFvNOsm0fZO2-ey2Z3mvtLjtoZzRbmLofGRaHO2F0qfdE14Bi1d-1hWY1dtlLYihIyw22ovrdKevEt1y2OomtC_CRg-JhSyTJJ49So0EzeSQ0cxc?AVOverride=1" -O libtorch.zip
RUN unzip libtorch.zip
WORKDIR /usr/share/SPOTLIGHT/spotlight_pt_port/
RUN ls -ltr --color=auto
RUN sed /LIBTORCH_LOC/d autobuild.sh > ab
RUN echo "LIBTORCH_LOC=/usr/share/libtorch" > autobuild.sh
RUN cat ab >> autobuild.sh
RUN chmod +x autobuild.sh && ./autobuild.sh
