#! /bin/sh

# get the input file
file=$1
echo ""
echo "You are checking file: $file"

# check if the input file is a Matlab script
if [[ "`echo $file | cut -d '.' -f2`" != "m" ]]; then
    echo "This file is not a Matlab script, quit."
    exit 1;
fi

###############
# Move the comment block below the function name
###############

# Find the comment block
count=1
flag=0
while read -r line
do
    first_c=`echo $line | cut -c 1`
    if [[ "$first_c" == "%" ]]; then
        if [[ $flag == 0 ]]; then
        start_num=$count
        flag=1
        # echo "start number: $start_num"
        fi
    fi
    if [[ $flag == 1 ]]; then
        if [[ "$first_c" != "%" ]]; then
            flag=0
            end_num=$(($count-1))
            # echo "end number: $end_num"
            break
        fi
    fi
    count=$(($count+1))
done < $file
length=$(($end_num-$start_num+1))
# echo "length: $length"

if [[ "$length" -lt "0" ]]; then
    echo -e "\e[31mNo comment block is found, please check file: $file"
    echo "Please add a comment block to explain what this function does"
    exit 1;
else

    # find the first function claim line
    echo ""
    echo "==> Start moving the comment block "

    line_num_array=`sed -n "/^function/=" ${file}`
    line_num=`echo $line_num_array | cut -c1`
    correct_start_line=$(($line_num+2))

    if [[ $correct_start_line == $start_num ]]; then
        echo "  [PASSED]: comment block is in correct format."
    else   
        echo -e "\e[31m  [ERROR] $file does not have a correct comment block's format\e[0m"

        # cp the comment block to a txt file
        sed -n "$start_num,$end_num"p $file > comment.txt

        echo ""

        # show the original file 
        echo -e "====================== Original file ==================="
        sed -n "1,$(($end_num+2))"p $file

        # show the new file with the fixed comment block

        echo ""
        echo -e "===================== Proposed file (with the comment block fixed) ===================="
        rsync -az $file ${file}_backup
        sed -i "$start_num,$end_num"d ${file}_backup
        line_num_array=`sed -n "/^function/=" ${file}_backup`
        line_num=`echo $line_num_array | cut -c1`
        sed -i "${line_num}G " ${file}_backup    
        sed -i "$((${line_num}+1))r comment.txt" ${file}_backup
        sed -n "1,$(($length+4))"p ${file}_backup

        # ask the user whether he wants to move this comment block
        echo ""
        read -r -p "===> Do want to fix the comment block as proposed above?[y/n]" response </dev/tty 
        if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
            mv ${file}_backup ${file}
        else
            rm ${file}_backup
            echo "  [Skipped]"
        fi
        rm comment.txt
        echo ""
    fi

    ###############
    # check and add function claim
    ###############
    echo ""
    echo "=> Start checking the function claim"
    RED=`echo -e '\033[31m'`
    BLUE=`echo -e '\033[34m'`
    NORMAL=`echo -e '\033[0m'`

    function_claim=`sed -n "$line_num"p $file`
    claim=${function_claim:9}
    # if the claim include bracket: [], replace it with \[ \]
    claim=$(echo $claim | sed "s/\[/\\\[/g" | sed "s/\]/\\\]/g")
    claim_block=`sed -n "$(($line_num+2)),$(($line_num+$length+1))"p $file | sed -n "/$claim/"p` 
    if [[ "$claim_block" != "" ]]; then
        echo "  [PASSED]"
    else
        echo -e "\e[31m  [ERROR] $file has no function claim\e[0m"
        echo -e "=================== Original file =================="
        sed -n "1,$(($length+4+$line_num))"p $file

        echo ""
        echo -e "==================== Proposed file (with the function claim added) =================="
        sed "$(($line_num+1))a\
            $BLUE% ${claim}$NORMAL" $file | head -n $(($length+4+$line_num)) 

        echo ""
        read -r -p "===> Do you want to add function claim as proposed above?[y/n]" response </dev/tty
    
        # whether you want to add a license or not
        if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
            sed -i "$(($line_num+1))a\
                % ${claim}" $file
            length=$(($length+1))
        else
            echo "  [Skipped]"
        fi
    fi


    ###############
    # check and add the license
    ###############

    echo ""
    echo "=> Start checking for license"

    license_line=`grep "MIT license" $file`
    #echo $license_line
    mit_license="% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md"

    # check whether there is the mit license
    if [[ "$license_line" == *"CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md" && "$license_line" == *"Written by"* ]] ; then
        echo -e "  [PASSED]"
    else 

        # check whether there is other license
        copyright=`sed -n '/copyright/p' $file`
        if [[ "$copyright" == "" ]]; then
            echo -e "\e[31m  [ERROR] $file has no license.\e[0m"

            echo "m===================== Original file =========================="
            sed -n "1,$(($length+4))"p $file

            echo ""
            echo -e "====================== New file (with MIT license added) ==========================="
            sed "$(($length+2+$line_num-1))a\%\n\
                $BLUE${mit_license}$NORMAL\n" $file | head -n $(($length+4+$line_num)) 

            echo ""
            read -r -p "===> Do you want to add the MIT license?[y/n]" response </dev/tty
        
            # whether you want to add a license or not
            if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
                sed -i "$(($length+2+$line_num-1))a\%\n\
                    ${mit_license}\n" $file
            else
                echo "  [Skipped]"
            fi
        else 
            echo -e "  [PASSED] \e[31m$file has a license: $copyright.\e[0m"
        fi
    fi
fi
