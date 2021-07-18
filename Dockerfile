FROM continuumio/miniconda3

COPY deleat_env.txt /home
COPY deleat-v0.1 /home/deleat-v0.1
RUN conda create --name deleat-v0.1 --file /home/deleat_env.txt
RUN echo "alias deleat='python /home/deleat-v0.1/deleat.py'" >> /etc/profile

CMD ["/bin/bash", "--login"]