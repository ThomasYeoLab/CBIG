#!/bin/sh
#
# Use as
#  package.sh <sourcedir> <destinationdir>

SOURCE=$1
DESTINATION=$2

if [ -e ${SOURCE} ]; then
SOURCE=${HOME}/matlab/fieldtrip
fi

if [ -e ${DESTINATION} ]; then
DESTINATION=${HOME}/matlab/cifti
fi

echo Copying from ${SOURCE} to ${DESTINATION}

# update the version
REV=`cd ${SOURCE} && git rev-parse HEAD`
date -u                                                   >  ${DESTINATION}/VERSION
echo https://github.com/fieldtrip/fieldtrip/commits/$REV  >> ${DESTINATION}/VERSION

# main functions
cp ${SOURCE}/fileio/ft_read_cifti.m                      ${DESTINATION} 
cp ${SOURCE}/fileio/ft_write_cifti.m                     ${DESTINATION} 

# external dependency
cp -r ${SOURCE}/external/gifti/@gifti                    ${DESTINATION} 
cp -r ${SOURCE}/external/gifti/@xmltree                  ${DESTINATION} 


# first order dependencies, i.e. called by main functions
cp ${SOURCE}/fileio/ft_read_headshape.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_write_headshape.m                 ${DESTINATION}/private
cp ${SOURCE}/fileio/private/fixname.m                    ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_getopt.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_getopt.mex*               ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_hastoolbox.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_warp_apply.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/getdimord.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/inflate_file.m               ${DESTINATION}/private
cp ${SOURCE}/fileio/private/istrue.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ndgrid.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/pos2transform.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_nifti2_hdr.m            ${DESTINATION}/private
cp ${SOURCE}/fileio/private/tokenize.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_nifti2_hdr.m           ${DESTINATION}/private
cp ${SOURCE}/utilities/copyfields.m                      ${DESTINATION}/private
cp ${SOURCE}/utilities/removefields.m                    ${DESTINATION}/private

# second order dependencies, i.e. called by first order functions
cp ${SOURCE}/fileio/ft_filetype.m                        ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_read_header.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_read_headshape.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_read_mri.m                        ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_read_sens.m                       ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_read_vol.m                        ${DESTINATION}/private
cp ${SOURCE}/fileio/ft_write_headshape.m                 ${DESTINATION}/private
cp ${SOURCE}/fileio/private/fetch_url.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/find_outermost_boundary.m    ${DESTINATION}/private
cp ${SOURCE}/fileio/private/fixname.m                    ${DESTINATION}/private
cp ${SOURCE}/fileio/private/fixpos.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_convert_units.m           ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_datatype_sens.m           ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_getopt.*                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_hastoolbox.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_scalingfactor.m           ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_senstype.m                ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_voltype.m                 ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_warp_apply.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/getdimord.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/inflate_file.m               ${DESTINATION}/private
cp ${SOURCE}/fileio/private/istrue.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ndgrid.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/pos2transform.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_asa.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_besa_sfp.m              ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_bti_hs.m                ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_bv_srf.m                ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_caret_spec.m            ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_ctf_hc.m                ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_ctf_pos.m               ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_ctf_shape.m             ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_neuromag_hc.m           ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_nifti2_hdr.m            ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_off.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_ply.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_polhemus_fil.m          ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_stl.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_vtk.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_yokogawa_header.m       ${DESTINATION}/private
cp ${SOURCE}/fileio/private/read_yokogawa_header_new.m   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/refine.m                     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/surf_to_tetgen.m             ${DESTINATION}/private
cp ${SOURCE}/fileio/private/tokenize.m                   ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_nifti2_hdr.m           ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_off.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_ply.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_stl.m                  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/write_vtk.m                  ${DESTINATION}/private
cp ${SOURCE}/utilities/copyfields.m                      ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_datatype.m                     ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_getopt.*                       ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_hastoolbox.m                   ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_scalingfactor.m                ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_struct2double.m                ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_warning.m                      ${DESTINATION}/private
cp ${SOURCE}/utilities/ft_warp_apply.m                   ${DESTINATION}/private
cp ${SOURCE}/utilities/hasyokogawa.m                     ${DESTINATION}/private
cp ${SOURCE}/utilities/istrue.m                          ${DESTINATION}/private
cp ${SOURCE}/utilities/keepfields.m                      ${DESTINATION}/private
cp ${SOURCE}/utilities/private/fixname.m                 ${DESTINATION}/private
cp ${SOURCE}/utilities/private/individual2sn.m           ${DESTINATION}/private
cp ${SOURCE}/utilities/private/sn2individual.m           ${DESTINATION}/private
cp ${SOURCE}/utilities/removefields.m                    ${DESTINATION}/private
cp ${SOURCE}/utilities/renamefields.m                    ${DESTINATION}/private
cp ${SOURCE}/utilities/tokenize.m                        ${DESTINATION}/private

# some other files that I noticed to be missing
cp ${SOURCE}/fileio/private/filetype_check_extension.m  ${DESTINATION}/private
cp ${SOURCE}/fileio/private/filetype_check_header.m     ${DESTINATION}/private
cp ${SOURCE}/fileio/private/filetype_check_uri.m        ${DESTINATION}/private
cp ${SOURCE}/fileio/private/ft_estimate_units.m         ${DESTINATION}/private
cp ${SOURCE}/fileio/private/getdimord.m                 ${DESTINATION}/private
cp ${SOURCE}/fileio/private/getdimsiz.m                 ${DESTINATION}/private
