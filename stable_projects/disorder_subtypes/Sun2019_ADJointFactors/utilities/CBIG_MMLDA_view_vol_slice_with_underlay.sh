#!/bin/sh

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# This script will overlay $in_vol on the $underlay_vol with $colormap and concat different slices horizontally.
# The cropping parameters assume FreeView's "1 big 3 small" mode, window fully expanded to 2000x1000, which means you
# need to create a VNC window with resolution 2000x1000. 

in_vol=$1
underlay_vol=$2
mni_space=$3
colormap=$4
plane=$5
out_dir=$6
out_name=$7

mkdir -p ${out_dir}

if [ "$mni_space" == "MNI2mm" ]; then
    if [ "$plane" == "coronal" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 20 45 -ss ${out_dir}/${out_name}_cor1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 30 45 -ss ${out_dir}/${out_name}_cor2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 40 45 -ss ${out_dir}/${out_name}_cor3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 50 45 -ss ${out_dir}/${out_name}_cor4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 60 45 -ss ${out_dir}/${out_name}_cor5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 70 45 -ss ${out_dir}/${out_name}_cor6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 45 80 45 -ss ${out_dir}/${out_name}_cor7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_cor1.png -crop +0+150 +repage -crop +0-50 +repage -crop +420+0 +repage -crop -420+0 +repage -flop ${out_dir}/${out_name}_cor1_shaved.png 
        convert ${out_dir}/${out_name}_cor2.png -crop +0+150 +repage -crop +0-50 +repage -crop +375+0 +repage -crop -375+0 +repage -flop ${out_dir}/${out_name}_cor2_shaved.png 
        convert ${out_dir}/${out_name}_cor3.png -crop +0+150 +repage -crop +0-50 +repage -crop +350+0 +repage -crop -350+0 +repage -flop ${out_dir}/${out_name}_cor3_shaved.png 
        convert ${out_dir}/${out_name}_cor4.png -crop +0+150 +repage -crop +0-50 +repage -crop +335+0 +repage -crop -335+0 +repage -flop ${out_dir}/${out_name}_cor4_shaved.png 
        convert ${out_dir}/${out_name}_cor5.png -crop +0+150 +repage -crop +0-50 +repage -crop +345+0 +repage -crop -345+0 +repage -flop ${out_dir}/${out_name}_cor5_shaved.png 
        convert ${out_dir}/${out_name}_cor6.png -crop +0+150 +repage -crop +0-50 +repage -crop +375+0 +repage -crop -375+0 +repage -flop ${out_dir}/${out_name}_cor6_shaved.png 
        convert ${out_dir}/${out_name}_cor7.png -crop +0+150 +repage -crop +0-50 +repage -crop +390+0 +repage -crop -390+0 +repage -flop ${out_dir}/${out_name}_cor7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_cor7_shaved.png \
        ${out_dir}/${out_name}_cor6_shaved.png \
        ${out_dir}/${out_name}_cor5_shaved.png \
        ${out_dir}/${out_name}_cor4_shaved.png \
        ${out_dir}/${out_name}_cor3_shaved.png \
        ${out_dir}/${out_name}_cor2_shaved.png \
        ${out_dir}/${out_name}_cor1_shaved.png \
        ${out_dir}/${out_name}_cor_concat.png

    elif [ "$plane" == "sagittal" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 75 54 45 -ss ${out_dir}/${out_name}_sag1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 65 54 45 -ss ${out_dir}/${out_name}_sag2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 55 54 45 -ss ${out_dir}/${out_name}_sag3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 45 54 45 -ss ${out_dir}/${out_name}_sag4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 35 54 45 -ss ${out_dir}/${out_name}_sag5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 25 54 45 -ss ${out_dir}/${out_name}_sag6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 15 54 45 -ss ${out_dir}/${out_name}_sag7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_sag1.png -crop +0+155 +repage -crop +0-55 +repage -crop +400+0 +repage -crop -440+0 +repage ${out_dir}/${out_name}_sag1_shaved.png 
        convert ${out_dir}/${out_name}_sag2.png -crop +0+155 +repage -crop +0-55 +repage -crop +325+0 +repage -crop -315+0 +repage ${out_dir}/${out_name}_sag2_shaved.png 
        convert ${out_dir}/${out_name}_sag3.png -crop +0+155 +repage -crop +0-55 +repage -crop +290+0 +repage -crop -290+0 +repage ${out_dir}/${out_name}_sag3_shaved.png 
        convert ${out_dir}/${out_name}_sag4.png -crop +0+155 +repage -crop +0-55 +repage -crop +295+0 +repage -crop -295+0 +repage ${out_dir}/${out_name}_sag4_shaved.png 
        convert ${out_dir}/${out_name}_sag5.png -crop +0+155 +repage -crop +0-55 +repage -crop +290+0 +repage -crop -290+0 +repage ${out_dir}/${out_name}_sag5_shaved.png 
        convert ${out_dir}/${out_name}_sag6.png -crop +0+155 +repage -crop +0-55 +repage -crop +325+0 +repage -crop -310+0 +repage ${out_dir}/${out_name}_sag6_shaved.png 
        convert ${out_dir}/${out_name}_sag7.png -crop +0+155 +repage -crop +0-55 +repage -crop +405+0 +repage -crop -420+0 +repage ${out_dir}/${out_name}_sag7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_sag7_shaved.png \
        ${out_dir}/${out_name}_sag6_shaved.png \
        ${out_dir}/${out_name}_sag5_shaved.png \
        ${out_dir}/${out_name}_sag4_shaved.png \
        ${out_dir}/${out_name}_sag3_shaved.png \
        ${out_dir}/${out_name}_sag2_shaved.png \
        ${out_dir}/${out_name}_sag1_shaved.png \
        ${out_dir}/${out_name}_sag_concat.png

    elif [ "$plane" == "axial" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 10 -ss ${out_dir}/${out_name}_axi1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 20 -ss ${out_dir}/${out_name}_axi2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 30 -ss ${out_dir}/${out_name}_axi3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 40 -ss ${out_dir}/${out_name}_axi4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 50 -ss ${out_dir}/${out_name}_axi5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 60 -ss ${out_dir}/${out_name}_axi6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 45 54 70 -ss ${out_dir}/${out_name}_axi7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_axi1.png -crop +0+55 +repage -crop +0-55 +repage -crop +415+0 +repage -crop -407+0 +repage ${out_dir}/${out_name}_axi1_shaved.png 
        convert ${out_dir}/${out_name}_axi2.png -crop +0+55 +repage -crop +0-55 +repage -crop +360+0 +repage -crop -355+0 +repage ${out_dir}/${out_name}_axi2_shaved.png 
        convert ${out_dir}/${out_name}_axi3.png -crop +0+55 +repage -crop +0-55 +repage -crop +335+0 +repage -crop -342+0 +repage ${out_dir}/${out_name}_axi3_shaved.png 
        convert ${out_dir}/${out_name}_axi4.png -crop +0+55 +repage -crop +0-55 +repage -crop +335+0 +repage -crop -335+0 +repage ${out_dir}/${out_name}_axi4_shaved.png 
        convert ${out_dir}/${out_name}_axi5.png -crop +0+55 +repage -crop +0-55 +repage -crop +340+0 +repage -crop -345+0 +repage ${out_dir}/${out_name}_axi5_shaved.png 
        convert ${out_dir}/${out_name}_axi6.png -crop +0+55 +repage -crop +0-55 +repage -crop +355+0 +repage -crop -365+0 +repage ${out_dir}/${out_name}_axi6_shaved.png 
        convert ${out_dir}/${out_name}_axi7.png -crop +0+55 +repage -crop +0-55 +repage -crop +430+0 +repage -crop -435+0 +repage ${out_dir}/${out_name}_axi7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_axi7_shaved.png \
        ${out_dir}/${out_name}_axi6_shaved.png \
        ${out_dir}/${out_name}_axi5_shaved.png \
        ${out_dir}/${out_name}_axi4_shaved.png \
        ${out_dir}/${out_name}_axi3_shaved.png \
        ${out_dir}/${out_name}_axi2_shaved.png \
        ${out_dir}/${out_name}_axi1_shaved.png \
        ${out_dir}/${out_name}_axi_concat.png
    fi
elif [ "$mni_space" == "MNI1.5mm" ]; then
    if [ "$plane" == "coronal" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 27 60 -ss ${out_dir}/${out_name}_cor1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 40 60 -ss ${out_dir}/${out_name}_cor2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 53 60 -ss ${out_dir}/${out_name}_cor3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 67 60 -ss ${out_dir}/${out_name}_cor4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 80 60 -ss ${out_dir}/${out_name}_cor5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 93 60 -ss ${out_dir}/${out_name}_cor6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport coronal -slice 60 107 60 -ss ${out_dir}/${out_name}_cor7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_cor1.png -crop +0+150 +repage -crop +0-50 +repage -crop +420+0 +repage -crop -420+0 +repage -flop ${out_dir}/${out_name}_cor1_shaved.png 
        convert ${out_dir}/${out_name}_cor2.png -crop +0+150 +repage -crop +0-50 +repage -crop +375+0 +repage -crop -375+0 +repage -flop ${out_dir}/${out_name}_cor2_shaved.png 
        convert ${out_dir}/${out_name}_cor3.png -crop +0+150 +repage -crop +0-50 +repage -crop +350+0 +repage -crop -350+0 +repage -flop ${out_dir}/${out_name}_cor3_shaved.png 
        convert ${out_dir}/${out_name}_cor4.png -crop +0+150 +repage -crop +0-50 +repage -crop +335+0 +repage -crop -335+0 +repage -flop ${out_dir}/${out_name}_cor4_shaved.png 
        convert ${out_dir}/${out_name}_cor5.png -crop +0+150 +repage -crop +0-50 +repage -crop +345+0 +repage -crop -345+0 +repage -flop ${out_dir}/${out_name}_cor5_shaved.png 
        convert ${out_dir}/${out_name}_cor6.png -crop +0+150 +repage -crop +0-50 +repage -crop +375+0 +repage -crop -375+0 +repage -flop ${out_dir}/${out_name}_cor6_shaved.png 
        convert ${out_dir}/${out_name}_cor7.png -crop +0+150 +repage -crop +0-50 +repage -crop +390+0 +repage -crop -390+0 +repage -flop ${out_dir}/${out_name}_cor7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_cor7_shaved.png \
        ${out_dir}/${out_name}_cor6_shaved.png \
        ${out_dir}/${out_name}_cor5_shaved.png \
        ${out_dir}/${out_name}_cor4_shaved.png \
        ${out_dir}/${out_name}_cor3_shaved.png \
        ${out_dir}/${out_name}_cor2_shaved.png \
        ${out_dir}/${out_name}_cor1_shaved.png \
        ${out_dir}/${out_name}_cor_concat.png

    elif [ "$plane" == "sagittal" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 100 72 60 -ss ${out_dir}/${out_name}_sag1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 87 72 60 -ss ${out_dir}/${out_name}_sag2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 73 72 60 -ss ${out_dir}/${out_name}_sag3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 60 72 60 -ss ${out_dir}/${out_name}_sag4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 47 72 60 -ss ${out_dir}/${out_name}_sag5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 33 72 60 -ss ${out_dir}/${out_name}_sag6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport sagittal -slice 20 72 60 -ss ${out_dir}/${out_name}_sag7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_sag1.png -crop +0+155 +repage -crop +0-55 +repage -crop +400+0 +repage -crop -440+0 +repage ${out_dir}/${out_name}_sag1_shaved.png 
        convert ${out_dir}/${out_name}_sag2.png -crop +0+155 +repage -crop +0-55 +repage -crop +325+0 +repage -crop -315+0 +repage ${out_dir}/${out_name}_sag2_shaved.png 
        convert ${out_dir}/${out_name}_sag3.png -crop +0+155 +repage -crop +0-55 +repage -crop +290+0 +repage -crop -290+0 +repage ${out_dir}/${out_name}_sag3_shaved.png 
        convert ${out_dir}/${out_name}_sag4.png -crop +0+155 +repage -crop +0-55 +repage -crop +295+0 +repage -crop -295+0 +repage ${out_dir}/${out_name}_sag4_shaved.png 
        convert ${out_dir}/${out_name}_sag5.png -crop +0+155 +repage -crop +0-55 +repage -crop +290+0 +repage -crop -290+0 +repage ${out_dir}/${out_name}_sag5_shaved.png 
        convert ${out_dir}/${out_name}_sag6.png -crop +0+155 +repage -crop +0-55 +repage -crop +325+0 +repage -crop -310+0 +repage ${out_dir}/${out_name}_sag6_shaved.png 
        convert ${out_dir}/${out_name}_sag7.png -crop +0+155 +repage -crop +0-55 +repage -crop +405+0 +repage -crop -420+0 +repage ${out_dir}/${out_name}_sag7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_sag7_shaved.png \
        ${out_dir}/${out_name}_sag6_shaved.png \
        ${out_dir}/${out_name}_sag5_shaved.png \
        ${out_dir}/${out_name}_sag4_shaved.png \
        ${out_dir}/${out_name}_sag3_shaved.png \
        ${out_dir}/${out_name}_sag2_shaved.png \
        ${out_dir}/${out_name}_sag1_shaved.png \
        ${out_dir}/${out_name}_sag_concat.png

    elif [ "$plane" == "axial" ]; then
        # overlay input vol on the underlay vol
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 13 -ss ${out_dir}/${out_name}_axi1
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 27 -ss ${out_dir}/${out_name}_axi2
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 40 -ss ${out_dir}/${out_name}_axi3
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 53 -ss ${out_dir}/${out_name}_axi4
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 67 -ss ${out_dir}/${out_name}_axi5
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 80 -ss ${out_dir}/${out_name}_axi6
        freeview -v ${underlay_vol} ${in_vol}:colormap=${colormap} -viewport axial -slice 60 72 93 -ss ${out_dir}/${out_name}_axi7

        # trim the black space of screenshorts
        convert ${out_dir}/${out_name}_axi1.png -crop +0+55 +repage -crop +0-55 +repage -crop +415+0 +repage -crop -407+0 +repage ${out_dir}/${out_name}_axi1_shaved.png 
        convert ${out_dir}/${out_name}_axi2.png -crop +0+55 +repage -crop +0-55 +repage -crop +360+0 +repage -crop -355+0 +repage ${out_dir}/${out_name}_axi2_shaved.png 
        convert ${out_dir}/${out_name}_axi3.png -crop +0+55 +repage -crop +0-55 +repage -crop +335+0 +repage -crop -342+0 +repage ${out_dir}/${out_name}_axi3_shaved.png 
        convert ${out_dir}/${out_name}_axi4.png -crop +0+55 +repage -crop +0-55 +repage -crop +335+0 +repage -crop -335+0 +repage ${out_dir}/${out_name}_axi4_shaved.png 
        convert ${out_dir}/${out_name}_axi5.png -crop +0+55 +repage -crop +0-55 +repage -crop +340+0 +repage -crop -345+0 +repage ${out_dir}/${out_name}_axi5_shaved.png 
        convert ${out_dir}/${out_name}_axi6.png -crop +0+55 +repage -crop +0-55 +repage -crop +355+0 +repage -crop -365+0 +repage ${out_dir}/${out_name}_axi6_shaved.png 
        convert ${out_dir}/${out_name}_axi7.png -crop +0+55 +repage -crop +0-55 +repage -crop +430+0 +repage -crop -435+0 +repage ${out_dir}/${out_name}_axi7_shaved.png

        # merge screenshorts of different slices horizontally
        convert -colorspace CMY +append \
        ${out_dir}/${out_name}_axi7_shaved.png \
        ${out_dir}/${out_name}_axi6_shaved.png \
        ${out_dir}/${out_name}_axi5_shaved.png \
        ${out_dir}/${out_name}_axi4_shaved.png \
        ${out_dir}/${out_name}_axi3_shaved.png \
        ${out_dir}/${out_name}_axi2_shaved.png \
        ${out_dir}/${out_name}_axi1_shaved.png \
        ${out_dir}/${out_name}_axi_concat.png
    fi
fi
