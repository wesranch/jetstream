

sudo docker container run -d -p 127.0.0.1:8000:7681 -v /home/wrancher/landis/input:/container/scripts esiil/landis_v8:latest


#update python

sudo docker cp run-landis.py zen_chaum:/container/scripts/

#or

sudo /usr/bin/docker exec zen_chaum: python3 container/scripts/run-landis.py