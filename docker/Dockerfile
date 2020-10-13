FROM nvidia/cuda:10.1-devel as ceo-build

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    pkg-config build-essential libssl-dev\
    noweb libcurl4-openssl-dev ca-certificates \
    curl llvm-dev libclang-dev clang

WORKDIR /usr/src/gicsdom/ceo/CEO
ADD ceo.tar.gz .
RUN make all

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /usr/src/gicsdom
ADD gicsdom.tar.gz .
RUN cargo build --release

FROM nvidia/cuda:10.1-runtime

ARG APP=/usr/src/gicsdom

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev 

ENV APP_USER=appuser

RUN groupadd $APP_USER \
    && useradd -g $APP_USER $APP_USER \
    && mkdir -p ${APP}

COPY --from=ceo-build /usr/src/gicsdom/target/release/fiducial /usr/local/bin/fiducial
COPY --from=ceo-build /usr/src/gicsdom/target/release/dome_seeing /usr/local/bin/dome_seeing
COPY --from=ceo-build /usr/src/gicsdom/target/release/glao_dome_seeing /usr/local/bin/glao_dome_seeing
COPY --from=ceo-build /usr/src/gicsdom/target/release/fiducial_atmosphere /usr/local/bin/fiducial_atmosphere

RUN chown -R $APP_USER:$APP_USER ${APP}

USER $APP_USER
WORKDIR ${APP}

CMD ["fiducial_atmosphere"]
