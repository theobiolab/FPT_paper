FROM ubuntu:20.04

ENV PYTHONPATH=/home/fpt/bin

WORKDIR /home/ftp

COPY . .

RUN chmod +x ftp_singularity_container.sh
RUN ./ftp_singularity_container.sh
RUN rm ftp_singularity_container.sh

EXPOSE 8888

CMD [ "jupyter", "lab", "--allow-root", "--no-browser", "--ip=0.0.0.0" ]