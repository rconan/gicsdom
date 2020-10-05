ceo: 
	cd ../CEO; git archive --format=tar.gz HEAD > ../docker/ceo.tar.gz

gicsdom:
	cd ../; git archive --format=tar.gz HEAD > docker/gicsdom.tar.gz

dock: Dockerfile
	docker build -t gmto.im/gicsdom .
