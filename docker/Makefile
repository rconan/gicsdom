AWS=docker run --rm -it -v ~/.aws:/root/.aws amazon/aws-cli

ceo:
	cd ../CEO; git archive --format=tar.gz HEAD > ../docker/ceo.tar.gz

gicsdom:
	cd ../; git archive --format=tar.gz HEAD > docker/gicsdom.tar.gz

dock: Dockerfile
	docker build -t gmto.im/gicsdom .

run:
	docker run --runtime=nvidia --rm --env-file=env_gicsdom gmto.im/gicsdom

push:
	$(AWS) ecr get-login-password --region us-west-2 | docker login --username AWS --password-stdin 378722409401.dkr.ecr.us-west-2.amazonaws.com
	docker tag gmto.im/gicsdom:latest 378722409401.dkr.ecr.us-west-2.amazonaws.com/gmto.im/gicsdom:latest
	docker push 378722409401.dkr.ecr.us-west-2.amazonaws.com/gmto.im/gicsdom:latest
