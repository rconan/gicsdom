FROM nvidia/cuda:10.1-devel as ceo-build

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    pkg-config build-essential libssl-dev\
    noweb libcurl4-openssl-dev \
    git wget bzip2 ca-certificates \
    curl llvm-dev libclang-dev clang

WORKDIR /usr/src/gicsdom/ceo/CEO
ADD ceo.tar.gz .
RUN make all

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /usr/src/gicsdom
ADD gicsdom.tar.gz .
RUN cargo build --release

#COPY ceo ceo/CEO/gmtMirrors
#COPY data .

#CMD ["./target/release/dome_seeing"]

FROM nvidia/cuda:10.1-runtime

ARG APP=/usr/src/gicsdom

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
            libcurl4-openssl-dev 

ENV APP_USER=appuser

RUN groupadd $APP_USER \
    && useradd -g $APP_USER $APP_USER \
    && mkdir -p ${APP}

COPY --from=ceo-build /usr/src/gicsdom/target/release/fiducial ${APP}/fiducial
COPY --from=ceo-build /usr/src/gicsdom/target/release/dome_seeing ${APP}/dome_seeing
COPY --from=ceo-build /usr/src/gicsdom/target/release/glao_dome_seeing ${APP}/glao_dome_seeing

RUN chown -R $APP_USER:$APP_USER ${APP}

USER $APP_USER
WORKDIR ${APP}

COPY ceo ceo/CEO/gmtMirrors
COPY data .
RUN mkdir cfd && mv CFD_CASES.yaml cfd/

CMD ["./glao_dome_seeing"]
