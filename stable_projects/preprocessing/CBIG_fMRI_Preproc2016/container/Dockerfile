# * Use multi-stage build to minimize the size of the final image

##################################################
# Stage 1: Install packages and copy files
##################################################
FROM centos:centos7.9.2009 AS tmp

# Install dependencies
RUN sed -i s/mirror.centos.org/vault.centos.org/g /etc/yum.repos.d/*.repo && \
    sed -i s/^#.*baseurl=http/baseurl=http/g /etc/yum.repos.d/*.repo && \
    sed -i s/^mirrorlist=http/#mirrorlist=http/g /etc/yum.repos.d/*.repo && \
    yum install -y unzip csh tcsh wget curl rsync which ca-certificates \
    libX11 libXcomposite libXcursor libXdamage libXext libXfixes libXft \
    libXi libXrandr libXrender libXScrnSaver libXt libXtst libXxf86vm \
    alsa-lib atk cairo cups-libs fontconfig GConf2 gtk2 gdk-pixbuf2 \
    gnome-vfs2 gstreamer1-plugins-base gstreamer1 pango libsndfile libxcb \
    libxslt git cmake make zlib-devel zlib-static libXp libXpm epel-release \
    openmotif gsl xorg-x11-fonts-misc PyQt4 netpbm-progs gnome-tweak-tool ed \
    libpng12 xorg-x11-server-Xvfb firefox libgfortran-static vim-common \
    lapack-devel lapack-static blas-devel blas-static python-devel \
    python3-devel libX11-devel libXmu-devel mesa-libGL-devel && \
    yum -y clean all

# Copy setup files
COPY stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/setup/ /app/setup/

# Install Conda environment and set up CBIG environment
RUN export PIP_DEFAULT_TIMEOUT=100 && bash /app/setup/CBIG_preproc_python_env_container_setup.sh && \
    echo 'source /app/setup/CBIG_config_in_container.sh' >> ~/.bashrc

# Copy other files
ENV CBIG_CODE_DIR=/app/CBIG_private

ARG SUB_DIR_TO_COPY=data/templates/volume/FSL_MNI152_masks
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

ARG SUB_DIR_TO_COPY=external_packages
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

ARG SUB_DIR_TO_COPY=setup
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

ARG SUB_DIR_TO_COPY=stable_projects/brain_parcellation/Yan2023_homotopic/parcellations/FreeSurfer
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

ARG SUB_DIR_TO_COPY=stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

COPY unit_tests/replace_unittest_flag ${CBIG_CODE_DIR}/unit_tests/

ARG SUB_DIR_TO_COPY=utilities
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

ARG SUB_DIR_TO_COPY=stable_projects/preprocessing/CBIG_fMRI_Preproc2016
COPY ${SUB_DIR_TO_COPY} ${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}

# Move files to the right place. Empty directories are created so that
# they can be properly mounted when running singularity containers)
ARG PROJECT_DIR=${CBIG_CODE_DIR}/${SUB_DIR_TO_COPY}
RUN mv ${PROJECT_DIR}/container/apps / && \
    mv ${PROJECT_DIR}/container/placeholder_git ${CBIG_CODE_DIR}/.git && \
    mkdir -p /apps/matlab/R2018b && mkdir -p /apps/matlab/MCR_R2018b/v95 && \
    mkdir -p ${PROJECT_DIR}/unit_tests/output && \
    mkdir -p /mnt/isilon/CSC1/Yeolab/CodeMaintenance/UnitTestData

# Modify exisiting scripts to fit the need of containerized preprocessing pipeline
RUN export PATH="/app/miniconda/bin:$PATH" && \
    source activate CBIG_py3 && \
    cp ${CBIG_CODE_DIR}/setup/CBIG_sample_config.sh /app/setup/CBIG_config_in_container.sh && \
    python ${PROJECT_DIR}/container/modify_existing.py

# Purge caches
RUN source /app/setup/CBIG_config_in_container.sh && conda activate CBIG_py3 && \
    conda clean -afy && pip cache purge

# Create entrypoint script
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'source /app/setup/CBIG_config_in_container.sh' >> /app/entrypoint.sh && \
    echo 'exec $PREPROC2016 "$@"' >> /app/entrypoint.sh

# Change permissions to 777 to ensure non-root user has root-like permission in singularity containers
RUN mkdir /input && mkdir /out
RUN chmod -R 777 /app /apps
RUN chmod -R 777 /home /out

##################################################
# Stage 2: Final stage to produce the image
##################################################
FROM centos:centos7.9.2009

COPY --from=tmp / /

##################################################
# Set ENTRYPOINT, CMD, and WORKDIR
##################################################
LABEL maintainer="Tian Fang tianfang@nus.edu.sg"
ENV CBIG_CODE_DIR=/app/CBIG_private
WORKDIR ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["-help"]
