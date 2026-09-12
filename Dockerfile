FROM kivy/buildozer:latest

USER root

RUN apt-get update && \
    apt-get install -y autopoint gettext libtool pkg-config automake autoconf
